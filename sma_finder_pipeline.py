"""Hail Batch (https://hail.is/docs/batch/service.html) pipeline that runs sma_finder.py in parallel on many samples."""

import hashlib
import os
import pandas as pd
import subprocess
import sys

from sma_finder import SMN_C840_POSITION_1BASED
from step_pipeline import pipeline, Backend, Localize, Delocalize, all_outputs_exist

DOCKER_IMAGE = "weisburd/sma_finder@sha256:ecec8bb6ec5b2d0027599f5677fa93dc8acd3b4cb67d78a78237b1f911ea9688"

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
    group.add_argument("--impute-genome-version-if-missing", action="store_true",
        help="For samples where the genome version isn't specified in the genome_version column, attempt to parse it "
             "out of the bam or cram header prior to running SMA Finder.")
    group.add_argument("--localize-via-copy", action="store_true",
        help="Localize read data via GSUTIL COPY rather than via CLOUDFUSE")
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

    # filter table to rows that have CRAM or BAM paths
    length_before = len(df)
    df = df[~df[args.cram_or_bam_path_column].isna() & ~df[args.crai_or_bai_path_column].isna()]
    num_filtered = length_before - len(df)
    if num_filtered > 0:
        print(f"Filtered out {num_filtered} out of {length_before} ({100*num_filtered/length_before:0.1f}%) of rows "
              f"where the cram or bam column or the crai or bai column were empty")

    # drop duplicate rows
    length_before = len(df)
    df = df.drop_duplicates(subset=[args.cram_or_bam_path_column, args.crai_or_bai_path_column])
    num_filtered = length_before - len(df)
    if num_filtered > 0:
        print(f"Filtered out {num_filtered} out of {length_before} ({100*num_filtered/length_before:0.1f}%) rows "
              f"because they had the same cram or bam path as previous rows (ie. were duplicates)")
    df = df.sort_values(args.sample_id_column)

    # validate genome_version_column arg if specified
    if args.genome_version_column:
        if args.genome_version_column in df.columns:
            invalid_genome_versions = set(df[args.genome_version_column]) - VALID_GENOME_VERSIONS
            if invalid_genome_versions:
                bad_genome_version_count = sum(~df[args.genome_version_column].isin(VALID_GENOME_VERSIONS))
                print(f"WARNING: The '{args.genome_version_column}' column in {args.sample_table} contains unexpected "
                      f"values: {invalid_genome_versions} in {bad_genome_version_count} out of {len(df)} rows. "
                      f"The only allowed values are: " + ", ".join(sorted(VALID_GENOME_VERSIONS)) + ". The unexpected "
                      "values will be cleared and those samples will be tested for both hg37 and hg38 coordinates")
                #df.loc[~df[args.genome_version_column].isin(VALID_GENOME_VERSIONS) : args.genome_version_column] = ""
        else:
            print(f"WARNING: {args.genome_version_column} column not found in {args.sample_table}. "
                  f"Will test each sample for both hg37 and hg38 coordinates.")

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

    original_df = pd.read_table(args.sample_table, dtype=str)

    steps = []
    print(f"Processing {len(df)} samples")
    for idx, row in df.iterrows():
        row_sample_id = row[args.sample_id_column]
        row_cram_or_bam_path = row[args.cram_or_bam_path_column]
        row_crai_or_bai_path = row[args.crai_or_bai_path_column]

        # decide which genome version(s) to use for this sample
        row_genome_version = row.get(args.genome_version_column)
        if row_genome_version not in VALID_GENOME_VERSIONS and args.impute_genome_version_if_missing:
            if args.verbose:
                print(f"Imputing genome version for {row_cram_or_bam_path}")
            row_genome_version = get_genome_version_from_bam_or_cram_header(row_cram_or_bam_path, args.gcloud_project)
            if not row_genome_version:
                print(f"WARNING: unable to impute genome version for {row_cram_or_bam_path}. Skipping...")
                continue

            original_df.loc[idx, args.genome_version_column] = row_genome_version
            original_df.to_csv(args.sample_table, index=False, header=True, sep="\t")
            if args.verbose:
                print(f"Saved {args.sample_table} with the imputed {row_genome_version} genome version for {row_cram_or_bam_path}")

        if row_genome_version not in VALID_GENOME_VERSIONS:
            # try all the reference genome versions
            print(f"Genome version not known for {row_sample_id} so will try reference genome versions:", ", ".join(REFERENCE_FASTA_PATH.keys()))
            row_genome_versions = list(REFERENCE_FASTA_PATH.keys())
        else:
            row_genome_versions = [row_genome_version]

        # process this sample using the decided upon genome version(s)
        for row_genome_version in row_genome_versions:
            output_dir = os.path.join(args.output_dir, row_genome_version)

            row_sample_type = row.get(args.sample_type_column)
            if row_sample_type:
                output_dir = os.path.join(output_dir, row_sample_type)

            # step1: run sma_finder.py
            s1 = bp.new_step(
                f"SMA pipeline: {row_sample_id}",
                arg_suffix="step1",
                image=DOCKER_IMAGE,
                cpu=0.25,
                memory="standard",
                output_dir=output_dir,
                delocalize_by=Delocalize.GSUTIL_COPY if args.localize_via_copy else Delocalize.COPY,
            )

            if args.localize_via_copy:
                s1.storage("75Gi")
                s1.switch_gcloud_auth_to_user_account()

            s1.command("set -euxo pipefail")

            reference_fasta_input, reference_fasta_fai_input = s1.inputs(
                REFERENCE_FASTA_PATH[row_genome_version],
                REFERENCE_FASTA_FAI_PATH[row_genome_version],
                localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

            # process input files
            localize_by = Localize.GSUTIL_COPY if args.localize_via_copy else Localize.HAIL_BATCH_CLOUDFUSE

            cram_or_bam_input = s1.input(row_cram_or_bam_path, localize_by=localize_by)
            crai_or_bai_input = s1.input(row_crai_or_bai_path, localize_by=localize_by)

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

            if len(row_genome_versions) > 1:
                # allow sma_finder.py command to fail when there is more than one reference genome to try
                s1.command("set +eo pipefail")

            # run smn_finder.py
            if row_genome_version.lower() != "t2t":
                row_genome_version_label = f"hg{row_genome_version}"
            else:
                row_genome_version_label = row_genome_version

            output_tsv_name = f"{OUTPUT_FILENAME_PREFIX}.{row_sample_id}" #.{row_genome_version_label}"
            if row_sample_type:
                output_tsv_name += f".{row_sample_type}"
            output_tsv_name += ".tsv"

            s1.command(
                f"time python3 -u /sma_finder.py "
                f"--{row_genome_version_label}-reference-fasta {reference_fasta_input} "
                f"--output-tsv {output_tsv_name} "
                f"--verbose "
                f"{local_bam_path}"
            )
            s1.command("ls")

            # delocalize the output tsv
            s1.output(output_tsv_name)
            steps.append(s1)

    # step2: combine tables from step1 into a single table
    s2 = bp.new_step(
        f"Combine {len(steps)} tables",
        image=DOCKER_IMAGE,
        cpu=1,
        memory="standard",
        output_dir=output_dir,
        delocalize_by=Delocalize.COPY,
        arg_suffix="step2",
    )
    s2.command("set -euxo pipefail")

    combined_output_tsv_filename = f"combined_results.{len(df)}_samples.{analysis_id}.tsv"
    for i, step in enumerate(steps):
        #if not all_outputs_exist(step):
        #    print(f"WARNING: skipping {step} step since its output(s) are missing")
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
    os.system(f"gsutil -m cp {os.path.join(output_dir, combined_output_tsv_filename)} .")
    result_df = pd.read_table(combined_output_tsv_filename)
    result_df.loc[:, "sample_id_or_filename"] = result_df.sample_id.fillna(result_df.filename_prefix)
    result_df = result_df.drop("filename_prefix", axis=1)

    df = df.drop_duplicates(subset=[args.sample_id_column], keep="first")
    df = df.drop(["genome_version"], axis=1)
    df_with_metadata = pd.merge(result_df, df, how="left", left_on="sample_id_or_filename", right_on=args.sample_id_column)
    df_with_metadata.to_csv(combined_output_tsv_filename, sep="\t", header=True, index=False)
    print(f"Wrote {len(df_with_metadata)} rows to {combined_output_tsv_filename}")


def get_genome_version_from_bam_or_cram_header(bam_or_cram_path, gcloud_project=None):
    gcloud_project_arg = f"-u {gcloud_project}" if gcloud_project else ""

    # get genome version from file header
    output = subprocess.check_output(
        f"gsutil {gcloud_project_arg} cat %s | samtools view -H - | grep @SQ | head -n 3" % bam_or_cram_path, shell=True, encoding="UTF-8", stderr=subprocess.DEVNULL)
    genome_version = None
    if "AS:GRCh37" in output or "Homo_sapiens_assembly19.fasta" in output:
        genome_version = "37"
    elif "AS:GRCh38" in output or "Homo_sapiens_assembly38.fasta" in output:
        genome_version = "38"
    else:
        print(f"WARNING: Unable to parse genome version from header lines in {bam_or_cram_path}: {output}")

    return genome_version


if __name__ == "__main__":
    main()
