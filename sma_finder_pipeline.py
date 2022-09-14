"""Hail Batch (https://hail.is/docs/batch/service.html) pipeline that runs sma_finder.py in parallel on many samples."""

import hashlib
import os
import pandas as pd
import sys

from sma_finder import SMN_C840_POSITION_1BASED
from step_pipeline import pipeline, Backend, Localize, Delocalize, all_outputs_exist

DOCKER_IMAGE = "weisburd/sma_finder@sha256:b9a635d0adce90a13fdbcdec4336dc2f124d50fe7c63074719379c8eb99a01f4"

REFERENCE_FASTA_PATH = {
    "37": "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta",
    "38": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
    #"T2T": None,
}
REFERENCE_FASTA_FAI_PATH = {
    "37": "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai",
    "38": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
    #"T2T": None,
}

VALID_GENOME_VERSIONS = set(REFERENCE_FASTA_PATH.keys())
VALID_SAMPLE_TYPES = {"WES", "WGS", "RNA"}

SMN_REGION_PADDING = 1000  # base-pairs

OUTPUT_FILENAME_PREFIX = "sma_finder"


def parse_args(batch_pipeline):
    """Define and parse command-line args"""

    arg_parser = batch_pipeline.get_config_arg_parser()
    group = arg_parser.add_argument_group("sma_finder general settings")
    group.add_argument("-o", "--output-dir", required=True,
        help="Output directory where to copy the sma_finder output .tsv file(s).")
    group.add_argument("-s", "--sample-id", action="append",
        help="If specified, only this sample id will be processed from the input table (useful for testing).")
    group.add_argument("-n", "--num-samples-to-process", type=int,
        help="If specified, only this many samples will be processed from the input table (useful for testing).")
    group.add_argument("-g", "--genome-version", choices=sorted(VALID_GENOME_VERSIONS),
        help="If all samples have the same genome version and the input table doesn't have a 'genome_version' column, "
             "it can be specified using this arg instead.")
    group.add_argument("-t", "--sample-type", choices=sorted(VALID_SAMPLE_TYPES, reverse=True),
        help="If all samples have the same sample type and the input table doesn't have a 'sample_type' column, "
             "it can be specified using this arg instead.")

    group.add_argument("sample_table",
        help="Path of tab-delimited table containing sample ids along with their BAM or CRAM file paths "
                        "and other metadata. The table should contain at least the following columns: "
                        "'genome_version', "
                        "'sample_type', "
                        "'sample_id' or 'individual_id', "
                        "'cram_path' or 'bam_path' or 'reads', "
                        "'bai_path' or 'crai_path' or 'index'. "
                        "Here, 'genome_version' must be one of: " + ", ".join(sorted(VALID_GENOME_VERSIONS)) + ", and "
                        "'sample_type' must be one of: " + ", ".join(sorted(VALID_SAMPLE_TYPES, reverse=True)))

    group = arg_parser.add_argument_group("sma_finder input table columns")
    group.add("--genome-version-column", default="genome_version",
              help="Optionally specify the name of input table column that contains the genome version: " +
                   ", ".join(sorted(VALID_GENOME_VERSIONS)))
    group.add("--sample-type-column", default="sample_type",
              help="Optionally specify the name of input table column that contains the sample type: " +
                   ", ".join(sorted(VALID_SAMPLE_TYPES, reverse=True)))
    group.add("--sample-id-column",
              help="Optionally specify the name of input table column that contains the sample id")
    group.add("--cram-or-bam-path-column",
              help="Optionally specify the name of input table column that contains the CRAM or BAM path")
    group.add("--crai-or-bai-path-column",
              help="Optionally specify the name of input table column that contains the CRAI or BAI path")
    args = batch_pipeline.parse_known_args()

    return args


def parse_sample_table(batch_pipeline):
    """Parse and validate the sample input table which contains paths and metadata for samples to process.

    Return:
        2-tuple (pandas.DataFrame, args): The input table and command line args.
    """
    args = parse_args(batch_pipeline)

    df = pd.read_table(args.sample_table, dtype=str)

    # check table columns
    arg_parser = batch_pipeline.get_config_arg_parser()
    if not args.sample_id_column:
        if "individual_id" in df.columns:
            args.sample_id_column = "individual_id"
        elif "sample_id" in df.columns:
            args.sample_id_column = "sample_id"
        else:
            arg_parser.error(f"{args.sample_table} must have one of these columns: 'individual_id', 'sample_id'")
    elif args.sample_id_column not in df.columns:
        arg_parser.error(f"{args.sample_table} doesn't have a '{args.sample_id_column}' column")

    if not args.cram_or_bam_path_column:
        if "cram_path" in df.columns:
            args.cram_or_bam_path_column = "cram_path"
        elif "bam_path" in df.columns:
            args.cram_or_bam_path_column = "bam_path"
        elif "reads" in df.columns:
            args.cram_or_bam_path_column = "reads"
        else:
            arg_parser.error(f"I{args.sample_table} must have one of these columns: 'cram_path', 'bam_path', 'reads'")
    elif args.cram_or_bam_path_column not in df.columns:
        arg_parser.error(f"{args.sample_table} doesn't have a '{args.cram_or_bam_path_column}' column")

    if not args.crai_or_bai_path_column:
        if "crai_path" in df.columns:
            args.crai_or_bai_path_column = "crai_path"
        elif "bai_path" in df.columns:
            args.crai_or_bai_path_column = "bai_path"
        elif "index" in df.columns:
            args.crai_or_bai_path_column = "index"
        else:
            arg_parser.error(f"{args.sample_table} must have one of these columns: 'crai_path', 'bai_path', 'index'")
    elif args.crai_or_bai_path_column not in df.columns:
        arg_parser.error(f"{args.sample_table} doesn't have a '{args.crai_or_bai_path_column}' column")

    if args.genome_version_column not in df.columns:
        if args.genome_version:
            df.loc[:, args.genome_version_column] = args.genome_version
        else:
            arg_parser.error(f"{args.sample_table} does not have a '{args.genome_version_column}' column. If all "
                             f"samples are aligned to the same reference genome, you can specify it using the "
                             f"--genome-version arg.")

    if args.sample_type_column not in df.columns:
        if args.sample_type:
            df.loc[:, args.sample_type_column] = args.sample_type
        else:
            arg_parser.error(f"{args.sample_table} does not have a '{args.sample_type_column}' column. If all samples "
                             f"are of the same type, you can specify it using the --sample-type arg.")

    # filter table to rows that have CRAM or BAM paths
    df = df[~df[args.cram_or_bam_path_column].isna() & ~df[args.crai_or_bai_path_column].isna()]
    df = df.drop_duplicates(subset=[args.cram_or_bam_path_column, args.crai_or_bai_path_column])
    df = df.sort_values(args.sample_id_column)

    # validate genome_version and sample_type columns
    invalid_genome_versions = set(df[args.genome_version_column]) - VALID_GENOME_VERSIONS
    if invalid_genome_versions:
        print(f"WARNING: The '{args.genome_version_column}' column in {args.sample_table} has unexpected values in some rows: "
              f"{invalid_genome_versions}. The only valid values are: " + ", ".join(sorted(VALID_GENOME_VERSIONS)) + ". ")
        length_before = len(df)
        df = df[~df[args.genome_version_column].isin(VALID_GENOME_VERSIONS)]
        print(f"Filtered out {len(df)} out of {length_before} ({100*len(df)/length_before:0.1f}%) of rows due to an "
              f"invalid '{args.genome_version_column}'")

    invalid_sample_types = set(df[args.sample_type_column]) - VALID_SAMPLE_TYPES
    if invalid_sample_types:
        print(f"WARNING: The '{args.sample_type_column}' column in {args.sample_table} has unexpected values in some "
              f"rows: {invalid_sample_types}. The only valid values are: " + ", ".join(sorted(VALID_SAMPLE_TYPES, reverse=True)) + ". ")
        length_before = len(df)
        df = df[~df[args.sample_type_column].isin(VALID_SAMPLE_TYPES)]
        print(f"Filtered out {len(df)} out of {length_before} ({100*len(df)/length_before:0.1f}%) of rows due to an "
              f"invalid '{args.sample_type_column}'")

    # apply --sample-id and --num-samples-to-process args if they were specified
    if args.sample_id:
        df = df[df[args.sample_id_column].isin(args.sample_id)]
        print(f"Found {len(df)} out of {len(args.sample_id)} requested sample ids")
        if len(df) == 0:
            sys.exit(1)
        elif len(df) < len(args.sample_id):
            print(f"WARNING: Couldn't find sample ids: {set(args.sample_id) - set(df.sample_id)}")

    if args.num_samples_to_process:
        df = df.iloc[:args.num_samples_to_process]

    print(f"Parsed {len(df)} rows from {args.sample_table}")

    return df, args


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE)

    df, args = parse_sample_table(bp)

    bp.set_name(f"sma_finder: {len(df)} samples")

    # compute a hash of the sample ids being processed
    analysis_id = ", ".join(df[args.sample_id_column])
    analysis_id = hashlib.md5(analysis_id.encode('UTF-8')).hexdigest().upper()
    analysis_id = analysis_id[:10]  # shorten

    if not args.force:
        bp.precache_file_paths(os.path.join(args.output_dir, "**", f"{OUTPUT_FILENAME_PREFIX}*"))

    steps = []
    print(f"Processing {len(df)} samples")
    for _, row in df.iterrows():
        row_sample_id = row[args.sample_id_column]
        row_genome_version = row[args.genome_version_column]
        row_sample_type = row[args.sample_type_column]
        row_cram_or_bam_path = row[args.cram_or_bam_path_column]
        row_crai_or_bai_path = row[args.crai_or_bai_path_column]

        # step1: run sma_finder.py
        s1 = bp.new_step(
            f"SMA pipeline: {row_sample_id}",
            arg_suffix="step1",
            image=DOCKER_IMAGE,
            cpu=0.25,
            memory="standard",
            output_dir=os.path.join(args.output_dir, row_genome_version, row_sample_type),
            delocalize_by=Delocalize.COPY,
        )
        s1.switch_gcloud_auth_to_user_account()
        s1.command("set -euxo pipefail")

        reference_fasta_input, reference_fasta_fai_input = s1.inputs(
            REFERENCE_FASTA_PATH[row_genome_version],
            REFERENCE_FASTA_FAI_PATH[row_genome_version],
            localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

        # process input files
        cram_or_bam_input = s1.input(row_cram_or_bam_path, localize_by=Localize.HAIL_BATCH_CLOUDFUSE_VIA_TEMP_BUCKET)
        crai_or_bai_input = s1.input(row_crai_or_bai_path, localize_by=Localize.HAIL_BATCH_CLOUDFUSE_VIA_TEMP_BUCKET)

        s1.command(f"ls -lh {cram_or_bam_input}")
        s1.command("cd /io/")

        # create symlinks in the same directory to work around cases when they are in different directories in the cloud
        s1.command(f"ln -s {cram_or_bam_input} /{cram_or_bam_input.filename}")
        s1.command(f"ln -s {crai_or_bai_input} /{crai_or_bai_input.filename}")

        # extract the regions of interest into a local bam file to avoid random access network requests downstream
        smn_chrom, smn1_position, _, smn2_position, _ = SMN_C840_POSITION_1BASED[row_genome_version]
        smn_interval_start = min(smn1_position, smn2_position) - SMN_REGION_PADDING
        smn_interval_end = max(smn1_position, smn2_position) + SMN_REGION_PADDING
        smn_region = f"{smn_chrom}:{smn_interval_start}-{smn_interval_end}"

        local_bam_path = f"{row_sample_id}.bam"
        s1.command(f"samtools view -T {reference_fasta_input} -b /{cram_or_bam_input.filename} {smn_region} "
                   f" | samtools sort > {local_bam_path}")
        s1.command(f"samtools index {local_bam_path}")

        # run smn_finder.py
        output_tsv_name = f"{OUTPUT_FILENAME_PREFIX}.{row_sample_id}.{row_sample_type}.tsv"
        s1.command(
            f"time python3 -u /sma_finder.py "
            f"-R {reference_fasta_input} "
            f"-g {row_genome_version} "
            f"--output-tsv {output_tsv_name} "
            f"--verbose "
            f"{local_bam_path}"
        )
        s1.command("ls")

        # delocalize the output tsv
        s1.output(output_tsv_name, delocalize_by=Delocalize.COPY)
        steps.append(s1)

    # step2: combine tables from step1 into a single table
    s2 = bp.new_step(
        f"Combine {len(steps)} tables",
        image=DOCKER_IMAGE,
        cpu=1,
        memory="standard",
        output_dir=args.output_dir,
        delocalize_by=Delocalize.COPY,
        arg_suffix="step2",
    )
    s2.command("set -euxo pipefail")

    combined_output_tsv_filename = f"combined_results.{len(df)}_samples.{analysis_id}.tsv"
    for i, step in enumerate(steps):
        #if args.skip_step1 and not all_outputs_exist(step):
        #    print(f"WARNING: skipping {step}")
        #    continue
        s2.depends_on(step)
        tsv_input = s2.use_previous_step_outputs_as_inputs(step, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
        if i == 0:
            s2.command(f"head -n 1 {tsv_input} > {combined_output_tsv_filename}")
        s2.command(f"tail -n +2 {tsv_input} >> {combined_output_tsv_filename}")

    s2.command(f"gzip {combined_output_tsv_filename}")
    combined_output_tsv_filename = f"{combined_output_tsv_filename}.gz"

    s2.output(combined_output_tsv_filename, delocalize_by=Delocalize.COPY)

    bp.run()

    # download the output table from step2 and merge it with the input table given to this pipeline.
    os.system(f"gsutil -m cp {os.path.join(args.output_dir, combined_output_tsv_filename)} .")
    result_df = pd.read_table(combined_output_tsv_filename)
    result_df.loc[:, "sample_id_or_filename"] = result_df.sample_id.where(
        result_df.sample_id.isin(set(df[args.sample_id_column])), result_df.filename)

    df = df.drop_duplicates(subset=[args.sample_id_column], keep="first")
    df_with_metadata = pd.merge(result_df, df, how="left", left_on="sample_id_or_filename", right_on=args.sample_id_column)
    df_with_metadata.to_csv(combined_output_tsv_filename, sep="\t", header=True, index=False)
    print(f"Wrote {len(df_with_metadata)} rows to {combined_output_tsv_filename}")


if __name__ == "__main__":
    main()
