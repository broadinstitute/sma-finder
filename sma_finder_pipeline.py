"""Hail Batch (https://hail.is/docs/batch/service.html) pipeline that runs sma_finder.py in parallel on many samples."""

import hashlib
import os
import pandas as pd
import re
import subprocess
import sys

from sma_finder import SMN_C840_POSITION_1BASED
from step_pipeline import pipeline, Backend, Localize, Delocalize, files_exist

DOCKER_IMAGE = "weisburd/sma_finder@sha256:101b94dc99ab0b17d18cd62db7cceb63c3e00f2561ae5459d537a31dacec0ceb"

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

SMN_REGION_PADDING = 1000  # base-pairs

OUTPUT_FILENAME_PREFIX = "sma_finder"


def parse_args(batch_pipeline):
    """Define and parse command-line args"""

    arg_parser = batch_pipeline.get_config_arg_parser()
    group = arg_parser.add_argument_group("pipeline settings")
    group.add_argument("-o", "--output-dir", required=True,
        help="Cloud storage output directory where to copy individual sma_finder output .tsv file(s).")
    group.add_argument("-s", "--sample-id", action="append",
        help="If specified, only this sample id will be processed from the input table (useful for testing).")
    group.add_argument("-n", "--num-samples-to-process", type=int,
        help="If specified, only this many samples will be processed from the input table (useful for testing).")
    group.add_argument("--offset", type=int, help="If specified, will skip his many initial rows from the table.")
    group.add_argument("-g", "--genome-version", choices=sorted(VALID_GENOME_VERSIONS),
        help="If all samples have the same genome version and the input table doesn't have a 'genome_version' column, "
             "it can be specified using this arg instead.")
    group.add_argument("--impute-genome-version-if-missing", action="store_true",
        help="For samples where the genome version isn't specified in the genome_version column, attempt to parse it "
             "out of the bam or cram header prior to running SMA Finder.")
    group.add_argument("--localize-via", choices={"cloudfuse", "copy", "gatk-print-reads"}, default="cloudfuse",
        help="Localize read data via CLOUDFUSE, COPY, or gatk PrintReads. Different options may be more or less "
             "cost-effective in different contexts")
    group.add_argument("--samples-per-job", type=int, default=1, help="Number of samples to process per Hail Batch "
        "Job. This is useful for reducing overhead involved in initializing each Job and localizing reference data.")

    group.add_argument("sample_table",
        help="Path of tab-delimited table containing sample ids along with their BAM or CRAM file paths "
                        "and other metadata. The table should contain at least the following columns: "
                        "'genome_version', "
                        "'sample_id' or 'individual_id', "
                        "'cram_path' or 'bam_path' or 'reads', "
                        "'bai_path' or 'crai_path' or 'index'. "
                        "Here, 'genome_version' must be one of: " + ", ".join(sorted(VALID_GENOME_VERSIONS)))

    group = arg_parser.add_argument_group("sma_finder input table columns")
    group.add("--genome-version-column", default="genome_version",
              help="Optionally specify the name of input table column that contains the genome version: " +
                   ", ".join(sorted(VALID_GENOME_VERSIONS)))
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
        print(f"Filtered out {num_filtered:,d} out of {length_before:,d} ({100*num_filtered/length_before:0.1f}%) of rows "
              f"where the cram or bam column or the crai or bai column were empty")

    # drop duplicate rows
    length_before = len(df)
    df = df.drop_duplicates(subset=[args.cram_or_bam_path_column, args.crai_or_bai_path_column])
    num_filtered = length_before - len(df)
    if num_filtered > 0:
        print(f"Filtered out {num_filtered:,d} out of {length_before:,d} ({100*num_filtered/length_before:0.1f}%) rows "
              f"because they had the same cram or bam path as previous rows (ie. were duplicates)")
    df = df.sort_values(args.sample_id_column)

    # validate genome_version_column arg if specified
    if args.genome_version_column:
        if args.genome_version_column in df.columns:
            invalid_genome_versions = set(df[args.genome_version_column]) - VALID_GENOME_VERSIONS
            if invalid_genome_versions:
                bad_genome_version_count = sum(~df[args.genome_version_column].isin(VALID_GENOME_VERSIONS))
                print(f"WARNING: The '{args.genome_version_column}' column in {args.sample_table} contains unexpected "
                      f"values: {invalid_genome_versions} in {bad_genome_version_count:,d} out of {len(df):,d} rows. "
                      f"The only allowed values are: " + ", ".join(sorted(VALID_GENOME_VERSIONS)) + ". The unexpected "
                      "values will be cleared and those samples will be tested for both hg37 and hg38 coordinates")
                #df.loc[~df[args.genome_version_column].isin(VALID_GENOME_VERSIONS) : args.genome_version_column] = ""
        else:
            print(f"WARNING: {args.genome_version_column} column not found in {args.sample_table}. "
                  f"Will test each sample for both hg37 and hg38 coordinates.")

    # apply --sample-id and --samples-to-process args if they were specified
    if args.sample_id:
        df = df[df[args.sample_id_column].isin(args.sample_id)]
        print(f"Found {len(df):,d} out of {len(args.sample_id):,d} requested sample ids")
        if len(df) == 0:
            sys.exit(1)
        elif len(df) < len(args.sample_id):
            print(f"WARNING: Couldn't find sample ids: {set(args.sample_id) - set(df.sample_id)}")

    print(f"Parsed {len(df):,d} rows from {args.sample_table}")

    return df, args


def output_tsv_path_func(output_dir, cram_path_column, sample_id_column):
    """This function returns a function that takes a row and returns an output tsv path"""
    def compute_output_tsv_path(row):
        cram_bucket_directories = os.path.dirname(re.sub("^gs://", "", row[cram_path_column]))
        sample_id = row[sample_id_column]
        return os.path.join(output_dir, f"sma_finder_tsvs/{cram_bucket_directories}",
                            f"{OUTPUT_FILENAME_PREFIX}.{sample_id}.tsv")

    return compute_output_tsv_path


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE)
    bp.get_config_arg_parser().add_argument("--skip-step2", action="store_true")
    bp.get_config_arg_parser().add_argument("--force-step2", action="store_true")    

    df, args = parse_sample_table(bp)
    if len(df) == 0:
        bp.get_config_arg_parser().error(f"{args.sample_table} doesn't contain any rows")

    if args.offset:
        if args.offset >= len(df):
            bp.get_config_arg_parser().error(
                f"{args.sample_table} contains only {len(df):,d} rows, which is smaller than the --offset value: {args.offset}")
        df = df.iloc[args.offset:]

    if args.num_samples_to_process:
        df = df.iloc[:args.num_samples_to_process]

    bp.set_name(f"sma_finder: {len(df):,d} samples")

    # compute output tsv
    df.loc[:, "output_tsv"] = df.apply(
        output_tsv_path_func(args.output_dir, args.cram_or_bam_path_column, args.sample_id_column), axis=1)

    df_copy = df.copy()

    # check for existing outputs
    if not args.force:
        existing_output_tsv_paths = {
            s["path"] for s in bp.precache_file_paths(os.path.join(args.output_dir, "**", f"{OUTPUT_FILENAME_PREFIX}*tsv"))
        }
        print(f"Found {len(existing_output_tsv_paths):,d} existing output tsv files in {args.output_dir}")
        before = len(df)
        df = df[~df["output_tsv"].isin(existing_output_tsv_paths)]
        print(f"Will skip {before - len(df):,d} samples that already have output tsv files")
    else:
        existing_output_tsv_paths = set()


    # compute a hash of the sample ids being processed
    analysis_id = ", ".join(sorted(df[args.sample_id_column]))
    analysis_id = hashlib.md5(analysis_id.encode('UTF-8')).hexdigest().upper()
    analysis_id = analysis_id[:10]  # shorten

    if args.genome_version_column in df.columns:
        df = df.sort_values(args.genome_version_column)


    if args.samples_per_job > 1:
        print(f"Processing {len(df):,d} samples in {int(len(df)/args.samples_per_job + 0.99):,d} batches of "
              f"{args.samples_per_job} sample(s) per job")
    else:
        print(f"Processing {len(df):,d} samples")

    for batch_start_i in range(0, len(df), args.samples_per_job):

        # create a new step1
        s1 = bp.new_step(
            step_number=1,
            arg_suffix="step1",
            image=DOCKER_IMAGE,
            cpu=0.25,
            memory="standard",
            delocalize_by=Delocalize.GSUTIL_COPY,
        )

        s1.name = "sma_finder: "
        df_current_batch = df.iloc[batch_start_i: batch_start_i + args.samples_per_job]
        if len(df_current_batch) > 3:
            s1.name += f"{len(df_current_batch):,d} samples"
        else:
            s1.name += ", ".join(df_current_batch[args.sample_id_column])

        # decide how to localize inputs
        if args.localize_via == "cloudfuse":
            localize_read_data_by = Localize.HAIL_BATCH_CLOUDFUSE
        elif args.localize_via == "copy":
            s1.storage("75Gi")
            s1.command("cd /io/")
            s1.switch_gcloud_auth_to_user_account()
            localize_read_data_by = Localize.GSUTIL_COPY
        elif args.localize_via == "gatk-print-reads":
            s1.switch_gcloud_auth_to_user_account(debug=True)
            localize_read_data_by = None
        else:
            raise ValueError(f"Unexpected localize_via arg: {args.localize_via}")

        localized_reference_fasta = {}
        samples_not_skipped_in_batch = 0
        for idx, row in df_current_batch.iterrows():
            row_sample_id = row[args.sample_id_column]
            row_cram_or_bam_path = row[args.cram_or_bam_path_column]
            row_crai_or_bai_path = row[args.crai_or_bai_path_column]
            row_output_tsv_path = row['output_tsv']

            # decide which genome version(s) to use for this sample
            row_genome_version = row.get(args.genome_version_column)
            if row_genome_version not in VALID_GENOME_VERSIONS and args.impute_genome_version_if_missing:
                if args.verbose:
                    print(f"Imputing genome version for {row_cram_or_bam_path}")
                row_genome_version = get_genome_version_from_bam_or_cram_header(row_cram_or_bam_path, args.gcloud_project)
                if not row_genome_version:
                    print(f"WARNING: unable to impute genome version for {row_cram_or_bam_path}. Skipping...")
                    continue

                df_copy.loc[idx, args.genome_version_column] = row_genome_version
                df_copy.to_csv(args.sample_table, index=False, header=True, sep="\t")
                if args.verbose:
                    print(f"Saved {args.sample_table} with the imputed {row_genome_version} genome version for {row_cram_or_bam_path}")

            if row_genome_version not in VALID_GENOME_VERSIONS:
                # try all the reference genome versions
                print(f"Genome version not known for {row_sample_id} so will try reference genome versions:", ", ".join(REFERENCE_FASTA_PATH.keys()))
                row_genome_versions = list(REFERENCE_FASTA_PATH.keys())
            else:
                row_genome_versions = [row_genome_version]


            if len(row_genome_versions) > 1:
                # allow sma_finder.py command to fail when there is more than one reference genome to try
                s1.command("echo set -ux >> command.sh")
            else:
                s1.command("echo set -euxo pipefail >> command.sh")

            s1.command("chmod 777 command.sh")

            # process this sample using the decided upon genome version(s)
            for row_genome_version in row_genome_versions:

                # localize the reference fasta if it hasn't been already
                if row_genome_version not in localized_reference_fasta:
                    localized_reference_fasta[row_genome_version] = s1.inputs(
                        REFERENCE_FASTA_PATH[row_genome_version],
                        REFERENCE_FASTA_FAI_PATH[row_genome_version],
                        localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

                reference_fasta_input, reference_fasta_fai_input = localized_reference_fasta[row_genome_version]

                # extract the regions of interest into a local bam file to avoid random access network requests downstream
                smn_chrom, smn1_position, _, smn2_position, _ = SMN_C840_POSITION_1BASED[row_genome_version]
                smn_region1 = f"{smn_chrom}:{smn1_position - SMN_REGION_PADDING}-{smn1_position + SMN_REGION_PADDING}"
                smn_region2 = f"{smn_chrom}:{smn2_position - SMN_REGION_PADDING}-{smn2_position + SMN_REGION_PADDING}"

                local_bam_path = f"{row_sample_id}.bam"

                if args.localize_via == "gatk-print-reads":
                    s1.command(f"echo 'gatk --java-options '-Xmx1G' PrintReads "
                               f"-R {reference_fasta_input} "
                               f"-L {smn_region1} "
                               f"-L {smn_region2} "
                               f"-I {row_cram_or_bam_path} "
                               f"-O {local_bam_path} "
                               f"--gcs-project-for-requester-pays {args.gcloud_project}' >> command.sh")
                else:
                    cram_or_bam_input = s1.input(row_cram_or_bam_path, localize_by=localize_read_data_by)
                    crai_or_bai_input = s1.input(row_crai_or_bai_path, localize_by=localize_read_data_by)

                    # create symlinks in the same directory to work around cases when they are in different directories in the cloud
                    s1.command(f"ln -s {cram_or_bam_input} /{cram_or_bam_input.filename}")
                    #s1.command(f"samtools index {cram_or_bam_input.filename} ")
                    s1.command(f"ln -s {crai_or_bai_input} /{crai_or_bai_input.filename}")

                    s1.command(f"echo 'samtools view -T {reference_fasta_input} -b /{cram_or_bam_input.filename} {smn_region1} {smn_region2} | samtools sort > {local_bam_path}' >> command.sh")

                s1.command(f"echo 'samtools index {local_bam_path}' >> command.sh")

                # run smn_finder.py
                if row_genome_version.lower() != "t2t":
                    row_genome_version_label = f"hg{row_genome_version}"
                else:
                    row_genome_version_label = row_genome_version

                samples_not_skipped_in_batch += 1
                s1.command(
                    f"echo 'python3 -u /sma_finder.py "
                    f"--{row_genome_version_label}-reference-fasta {reference_fasta_input} "
                    f"--output-tsv {os.path.basename(row_output_tsv_path)} "
                    f"--verbose "
                    f"{local_bam_path}' >> command.sh"
                )

                s1.command("./command.sh || true")
                s1.command("rm command.sh")

                # delocalize the output tsv
                s1.output(os.path.basename(row_output_tsv_path), row_output_tsv_path)

                s1.command(f"rm {local_bam_path} {local_bam_path}.bai")
                if args.localize_via == "copy":
                    # delete the local bam file
                    s1.command(f"rm {cram_or_bam_input.local_path} {crai_or_bai_input.local_path}")

            if samples_not_skipped_in_batch == 0:
                s1.cancel()
    bp.run()

    existing_output_tsv_paths = {
        s["path"] for s in bp.precache_file_paths(os.path.join(args.output_dir, "**", f"{OUTPUT_FILENAME_PREFIX}*tsv"))
    }
    print(f"Found {len(existing_output_tsv_paths):,d} output tsv files in {args.output_dir}")

    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE)
    bp.get_config_arg_parser().add_argument("--skip-step1", action="store_true")
    bp.get_config_arg_parser().add_argument("--force-step1", action="store_true")    
    args2 = parse_args(bp)

    df = df_copy
    df = df[df["output_tsv"].isin(existing_output_tsv_paths)]

    bp.set_name(f"sma_finder: combine {len(df)} samples")
    
    # step2: combine tables from step1 into a single table
    s2 = bp.new_step(
        f"Combine {len(df):,d} tables",
        image=DOCKER_IMAGE,
        cpu=1,
        memory="standard",
        output_dir=args.output_dir,
        localize_by=Localize.HAIL_BATCH_CLOUDFUSE, #GSUTIL_COPY,
        delocalize_by=Delocalize.COPY,
        step_number=2,
        arg_suffix="step2",
    )
    s2.switch_gcloud_auth_to_user_account()    
    s2.command("set -euxo pipefail")

    combined_output_tsv_filename = f"combined_results.{len(df)}_samples.{analysis_id}.tsv"
    for i, input_tsv_path in enumerate(set(df["output_tsv"])):
        input_tsv = s2.input(input_tsv_path)
        if i == 0:
            s2.command(f"head -n 1 {input_tsv} > {combined_output_tsv_filename}")
        s2.command(f"tail -n +2 {input_tsv} >> {combined_output_tsv_filename}")

    s2.command(f"gzip {combined_output_tsv_filename}")
    combined_output_tsv_filename = f"{combined_output_tsv_filename}.gz"

    s2.output(combined_output_tsv_filename, delocalize_by=Delocalize.COPY)

    bp.run()

    # download the output table from step2 and merge it with the input table given to this pipeline.
    os.system(f"gsutil -m cp {os.path.join(args.output_dir, combined_output_tsv_filename)} .")
    result_df = pd.read_table(combined_output_tsv_filename)
    result_df.loc[:, "sample_id_or_filename"] = result_df[args.sample_id_column].fillna(result_df.filename_prefix)
    result_df = result_df.drop("filename_prefix", axis=1)

    #df = df.drop_duplicates(subset=[args.sample_id_column], keep="first")
    df = df.drop(["genome_version"], axis=1)
    df_with_metadata = pd.merge(result_df, df, how="left", left_on="sample_id_or_filename", right_on=args.sample_id_column)
    df_with_metadata.to_csv(combined_output_tsv_filename, sep="\t", header=True, index=False)
    print(f"Wrote {len(df_with_metadata):,d} rows to {combined_output_tsv_filename}")


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
