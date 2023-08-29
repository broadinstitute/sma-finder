"""This script checks SMA Finder results for new samples with a 'has_sma' call. It takes a combined results table
which is created by running SMA Finder on multiple samples and concatenating the results into a single table.
"""

import argparse
import os
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("-m", "--metadata-table", action="append", help="optionally specify a sample metadata table(s) to join "
               "with results")
p.add_argument("-c", "--metadata-sample-id-column", default="sample_id", help="Sample id column in the metadata table")
p.add_argument("-k", "--known-case-sample-ids-path", help="Path of a text file containing a list of sample ids (one "
               "per line) that are already known to have SMA")
p.add_argument("-x", "--exclude-sample-ids-path", help="Path of a text file containing a list of sample ids (one per "
               "line) to exclude from the results")
p.add_argument("results_table", help="Path of the combined results table created by running SMA Finder on multiple "
               "samples and concatenating the output tables.")
args = p.parse_args()

# parse the results table
if not os.path.isfile(args.results_table):
    p.error(f"File not found: {args.results_table}")

df = pd.read_table(args.results_table)
df_standard_columns = set(df.columns)
print(f"Parsed {len(df):,d} rows from {args.results_table}")

def process_sample_id(s):
    return s.replace("-", "_").replace(" ", "_").replace(".", "_")

df["processed_sample_id"] = df["filename_prefix"].apply(process_sample_id)
df.set_index("processed_sample_id", inplace=True)

# join with metadata table if one was specified
df_metadata = None
if args.metadata_table:
    for metadata_table in args.metadata_table:
        if not os.path.exists(metadata_table):
            p.error(f"File not found: {metadata_table}")

    for metadata_table in args.metadata_table:
        current_df_metadata = pd.read_table(metadata_table)
        print(f"Parsed {len(current_df_metadata):,d} rows from {metadata_table}")
        if args.metadata_sample_id_column not in current_df_metadata.columns:
            p.error(f"--metadata-sample-id-column '{args.metadata_sample_id_column}' not found in sample metadata table "
                    f"{metadata_table}")
        current_df_metadata[args.metadata_sample_id_column] = current_df_metadata[args.metadata_sample_id_column].apply(
            process_sample_id)
        current_df_metadata.set_index(args.metadata_sample_id_column, inplace=True)
        if df_metadata is None:
            df_metadata = current_df_metadata
        else:
            df_metadata = pd.concat([df_metadata, current_df_metadata])

    sample_ids_with_metadata = set(df.index) & set(df_metadata.index)
    print(f"Found metadata entries for {len(sample_ids_with_metadata):,d} out of {len(set(df.index)):,d} "
          f"({len(sample_ids_with_metadata)/len(set(df.index)):.1%}) sample ids in {args.results_table}")

    df = df.join(df_metadata, how="left", rsuffix=":meta")

# load known case sample ids
known_cases = set()
if args.known_case_sample_ids_path:
    if not os.path.exists(args.known_case_sample_ids_path):
        p.error(f"File not found: {args.known_case_sample_ids_path}")

    with open(args.known_case_sample_ids_path) as f:
        known_cases = {
            process_sample_id(line.strip()) for line in f if line.strip()
        }
    print(f"Loaded {len(known_cases):,d} known case sample ids from {args.known_case_sample_ids_path}")


# load excluded sample ids
exclude_sample_ids = set()
if args.exclude_sample_ids_path:
    if not os.path.exists(args.exclude_sample_ids_path):
        p.error(f"File not found: {args.exclude_sample_ids_path}")

    with open(args.exclude_sample_ids_path) as f:
        exclude_sample_ids_path = {
            process_sample_id(line.strip()) for line in f if line.strip()
        }
    print(f"Loaded {len(exclude_sample_ids_path):,d} known case sample ids from {args.exclude_sample_ids_path}")

    df = df[~df.index.isin(exclude_sample_ids_path)]

# print result stats
sample_ids_not_has_sma = set(df[df.sma_status != 'has SMA'].index)
sample_ids_has_sma = set(df[df.sma_status == 'has SMA'].index)

print(f"SMA finder correctly called {len(sample_ids_has_sma & known_cases):,d} and missed "
      f"{len(sample_ids_not_has_sma & known_cases):,d} out of the {len(known_cases):,d} known case sample ids that "
      f"were included in this callset")
known_cases_outside_the_callset = known_cases - set(df.index)
if known_cases_outside_the_callset:
    print(f"The following known cases were not included in this callset: {', '.join(sorted(known_cases_outside_the_callset))}")
print(f"Overall, {len(sample_ids_has_sma):,d} samples have SMA in this callset: {args.results_table}")

df = df[df.sma_status == "has SMA"]
df = df[~df.index.isin(known_cases)]

# print stats on new cases
if len(df) > 0:
    print(f"\nNew cases:")
    extra_columns_to_print = set(df.columns) - df_standard_columns
    extra_columns_to_print = {c for c in extra_columns_to_print if "path" not in c}
    extra_columns_to_print = list(sorted(extra_columns_to_print))

    column_widths = {c: len(c) for c in extra_columns_to_print}
    for i, (_, row) in enumerate(df.iterrows()):
        for c in extra_columns_to_print:
            column_widths[c] = max(column_widths.get(c, 0), len(str(row.get(c, ""))))
    header_line = f"{'#':<2s}  {'filename_prefix':20s}  {'genome_version':20s}  "
    header_line += f"read support" + " "*15
    for c in extra_columns_to_print:
        header_line += f"{c:>{column_widths[c] + 5}s}"

    print(header_line)
    for i, (_, row) in enumerate(df.iterrows()):
        print_line = f"#{i+1:<2d}  {row['filename_prefix']:20s}  {row['genome_version']:20s}  "
        print_line += f"{int(row['c840_reads_with_smn1_base_C']):5,d} out of {int(row['c840_total_reads']):5,d} reads"
        row_dict = row.to_dict()
        for c in extra_columns_to_print:
            print_line += f" {str(row_dict.get(c, '')):>{column_widths[c] + 5}s}"
        print(print_line)
else:
    print("No new cases")
