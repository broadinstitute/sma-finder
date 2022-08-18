#!/usr/env python3

"""This script computes read counts at key positions in the SMN1 & SMN2 genes and uses these to determine whether
a sample has spinal muscular atrophy (SMA).
"""

import argparse
import os
import pysam
import pandas as pd
import pprint
import re

SMN_CHROMOSOME = {
    "37": "5",
    "38": "chr5",
    "T2T": "5",
}

VALID_GENOME_VERSIONS = set(SMN_CHROMOSOME.keys())

SMN_DIFFERING_POSITION_1BASED = {
    "37": {
        "c840":   {"SMN1": 70247773, "SMN2": 69372353},
        "c1124":  {"SMN1": 70248501, "SMN2": 69373081},
    },
    "38": {
        "c840":   {"SMN1": 70951946, "SMN2": 70076526},
        "c1124":  {"SMN1": 70952674, "SMN2": 70077254},
    },
    "T2T": {
        "c840":   {"SMN1": 71408734, "SMN2": 70810812},
        "c1124":  {"SMN1": 71409463, "SMN2": 70810084},
    }
}

SMN_OTHER_EXON_POSITIONS_1BASED = {
    "37": {
        "exon1":  {"SMN1": 70220962, "SMN2": 69345544},
        "exon2a": {"SMN1": 70234701, "SMN2": 69359277},
        "exon2b": {"SMN1": 70237275, "SMN2": 69361851},
        "exon3":  {"SMN1": 70238285, "SMN2": 69362861},
        "exon4":  {"SMN1": 70238621, "SMN2": 69363197},
        "exon5":  {"SMN1": 70240532, "SMN2": 69365109},
        "exon6":  {"SMN1": 70241948, "SMN2": 69366523},
    },
    "38": {
        "exon1":  {"SMN1": 70925135, "SMN2": 70049717},
        "exon2a": {"SMN1": 70938874, "SMN2": 70063450},
        "exon2b": {"SMN1": 70941448, "SMN2": 70066024},
        "exon3":  {"SMN1": 70942458, "SMN2": 70067034},
        "exon4":  {"SMN1": 70942794, "SMN2": 70067370},
        "exon5":  {"SMN1": 70944705, "SMN2": 70069282},
        "exon6":  {"SMN1": 70946121, "SMN2": 70070696},
    },
    "T2T":{
        "exon1":  {"SMN1": 71381923, "SMN2": 70837627},
        "exon2a": {"SMN1": 71395666, "SMN2": 70823885},
        "exon2b": {"SMN1": 71398240, "SMN2": 70821311},
        "exon3":  {"SMN1": 71399250, "SMN2": 70820301},
        "exon4":  {"SMN1": 71399586, "SMN2": 70819965},
        "exon5":  {"SMN1": 71401496, "SMN2": 70818055},
        "exon6":  {"SMN1": 71402910, "SMN2": 70816641},
    },
}

"""This threshold is based on the empirical distribution in positive and negative control samples from large 
WES and WGS cohorts.
"""
MIN_COVERAGE_NEEDED_TO_CALL_SMA_STATUS = 15
MAX_FRACTION_OF_READS_WITH_ERROR = 0.02
MAX_READS_WITH_ERROR = 2


def parse_args():
    """Parse and return command-line args"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference genome FASTA file path")
    parser.add_argument("-g", "--genome-version", required=True, choices=sorted(VALID_GENOME_VERSIONS),
                        help="Reference genome version")
    parser.add_argument("-o", "--output-tsv", help="Output tsv file path", default="output.tsv")
    parser.add_argument("-v", "--verbose", action="store_true", help="Whether to print extra details during the run")
    parser.add_argument("cram_or_bam_path", nargs="+", help="One or more CRAM or BAM file paths")
    args = parser.parse_args()

    # make sure the input files exist
    for path in [args.reference_fasta] + args.cram_or_bam_path:
        if not os.path.isfile(path):
            parser.error(f"File not found: {path}")

    return args


def count_nucleotides_at_position(alignment_file, chrom, pos_1based):
    """Count the number of A, C, G, and T's found at the given genomic position within the given read data.

    Args:
        alignment_file (pysam.AlignmentFile): The pysam AlignmentFile object representing the input BAM or CRAM file.
        chrom (str): Chromosome of the genomic position.
        pos_1based (int): 1-based genomic position where to count nucleotides.

    Return:
        dict: The keys are nucleotides "A", "C", "G", "T", and the values are counters representing the number of times
            a read contained that nucleotide at the given position.
    """

    pos_0based = pos_1based - 1
    nucleotide_counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    for pileup_column in alignment_file.pileup(
            region=f"{chrom}:{pos_0based}-{pos_1based}",
            stepper="samtools",
            ignore_overlaps=True,
            ignore_orphans=False,
            min_mapping_quality=0,
            min_base_quality=13,
            truncate=True,
            multiple_iterators=False):
        if pileup_column.pos < pos_0based:
            continue

        if pileup_column.pos != pos_0based:
            raise ValueError(f"Unexpected pileup position: {chrom}:{pileup_column.pos}. "
                             f"Expecting {chrom}:{pos_0based}")

        for base in pileup_column.get_query_sequences():
            if not base:
                continue

            base = base.upper()
            if base in nucleotide_counts:
                nucleotide_counts[base] += 1
            else:
                raise ValueError(f"Unexpected base '{base}' found at {chrom}:{pos_1based}")
        break

    return nucleotide_counts


def add_filename_and_file_type(cram_or_bam_path, output_row):
    """Add 'filename' and 'file_type' fields to the output_row.

    Args:
        cram_or_bam_path (str): Input CRAM or BAM path.
        output_row (dict): fields that will be written to the output .tsv
    """
    filename_pattern_match = re.match("(.*)[.](cram|bam)$", os.path.basename(cram_or_bam_path))
    if not filename_pattern_match:
        raise ValueError(f"File path doesn't end with have '.cram' or '.bam': {cram_or_bam_path}")

    output_row.update({
        "filename": filename_pattern_match.group(1),
        "file_type": filename_pattern_match.group(2),
    })


def determine_sample_id(alignment_file, output_row):
    """Try to get the sample_id from the BAM file header by looking for a read group (@RG) with a sample (SM) field.

    Args:
        alignment_file (pysam.AlignmentFile): The pysam AlignmentFile object representing the input BAM or CRAM file.
        output_row (dict): fields that will be written to the output .tsv
    """
    output_row["sample_id"] = ""
    for read_group in alignment_file.header.get("RG", []):
        if "SM" in read_group:
            output_row["sample_id"] = read_group["SM"]
            break


def count_reads_at_differing_bases(alignment_file, genome_version, output_row):
    """Count reads at SMN_DIFFERING_POSITION_1BASED and 'c840_reads_with_sm1_base_C', 'c840_total_reads', ... fields
    to the output_row.

    Args:
        alignment_file (pysam.AlignmentFile): The pysam AlignmentFile object representing the input BAM or CRAM file.
        genome_version (str): sample genome version (eg. "38")
        output_row (dict): fields that will be written to the output .tsv
    """
    chrom = SMN_CHROMOSOME[genome_version]
    smn_differing_positions = SMN_DIFFERING_POSITION_1BASED[genome_version]
    for key, smn1_base in [("c840", "C"), ("c1124", "G")]:
        differing_position = smn_differing_positions[key]
        smn1_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, differing_position["SMN1"])
        smn2_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, differing_position["SMN2"])
        output_row.update({
            f"{key}_reads_with_smn1_base_{smn1_base}": smn1_nucleotide_counts[smn1_base] + smn2_nucleotide_counts[smn1_base],
            f"{key}_total_reads": sum(smn1_nucleotide_counts.values()) + sum(smn2_nucleotide_counts.values()),
        })


def count_reads_at_other_exons(alignment_file, genome_version, output_row):
    """Compute coverage at SMN_OTHER_EXON_POSITIONS_1BASED and add 'reads_at_exon1', 'reads_at_exon2a', ... fields
    to the output_row.

    Args:
        alignment_file (pysam.AlignmentFile): The pysam AlignmentFile object representing the input BAM or CRAM file.
        genome_version (str): sample genome version (eg. "38")
        output_row (dict): fields that will be written to the output .tsv
    """

    chrom = SMN_CHROMOSOME[genome_version]
    for exon_name, exon_positions in SMN_OTHER_EXON_POSITIONS_1BASED[genome_version].items():
        smn1_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, exon_positions["SMN1"])
        smn2_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, exon_positions["SMN2"])
        output_row.update({
            f"reads_at_{exon_name}": sum(smn1_nucleotide_counts.values()) + sum(smn2_nucleotide_counts.values()),
        })


def call_sma_status(genome_version, output_row):
    """Determines SMA status and then sets the 'sma_status', 'sma_status_details', and 'average_exon_coverage' fields in
    output_row.

    Args:
        genome_version (str): sample genome version (eg. "38")
        output_row (dict): fields that will be written to the output .tsv
    """

    exon_positions = SMN_OTHER_EXON_POSITIONS_1BASED[genome_version]
    average_coverage_of_other_exons = sum(output_row[f"reads_at_{exon_name}"] for exon_name in exon_positions)
    average_coverage_of_other_exons /= len(exon_positions)

    if output_row['c840_total_reads'] >= MIN_COVERAGE_NEEDED_TO_CALL_SMA_STATUS:
        if (
            output_row['c840_reads_with_smn1_base_C'] <= MAX_READS_WITH_ERROR and
            output_row['c840_reads_with_smn1_base_C'] <= MAX_FRACTION_OF_READS_WITH_ERROR * output_row['c840_total_reads']
        ):
            sma_status = "has SMA"
            sma_status_details = "has 0 copies of SMN1"
        else:
            sma_status = "doesn't have SMA"
            sma_status_details = "has 1 or more copies of SMN1"
    else:
        if (
            average_coverage_of_other_exons >= MIN_COVERAGE_NEEDED_TO_CALL_SMA_STATUS and
            output_row['c840_reads_with_smn1_base_C'] <= MAX_READS_WITH_ERROR and
            output_row['c840_reads_with_smn1_base_C'] <= MAX_FRACTION_OF_READS_WITH_ERROR * average_coverage_of_other_exons
        ):
            sma_status = "may have SMA"
            sma_status_details = (
                    f"SMN exons 1-6 have coverage >= {MIN_COVERAGE_NEEDED_TO_CALL_SMA_STATUS}x while exon 7 has "
                    f"approximately zero coverage, suggesting either technical problems with sequencing exon 7, or "
                    f"alternatively the deletion of exon 7 in all copies of SMN1 and SMN2"
            )

        else:
            sma_status = "not enough coverage"
            sma_status_details = f"average coverage of all SMN exons is < {MIN_COVERAGE_NEEDED_TO_CALL_SMA_STATUS}x"

    output_row["sma_status"] = sma_status
    output_row["sma_status_details"] = sma_status_details
    output_row["average_exon_coverage"] = average_coverage_of_other_exons


def main():
    args = parse_args()

    # process the input BAM or CRAM files
    output_rows = []
    for cram_or_bam_path in args.cram_or_bam_path:
        output_row = {}
        add_filename_and_file_type(cram_or_bam_path, output_row)

        alignment_file = pysam.AlignmentFile(cram_or_bam_path, 'rc', reference_filename=args.reference_fasta)
        determine_sample_id(alignment_file, output_row)
        count_reads_at_differing_bases(alignment_file, args.genome_version, output_row)
        count_reads_at_other_exons(alignment_file, args.genome_version, output_row)

        call_sma_status(args.genome_version, output_row)

        output_rows.append(output_row)
        if args.verbose:
            print("Output row:")
            pprint.pprint(output_row)

        alignment_file.close()

    # write results to .tsv
    df = pd.DataFrame(output_rows)
    reordered_output_columns = ["filename", "file_type", "sample_id", "sma_status", "sma_status_details"]
    reordered_output_columns += [c for c in df.columns if c not in reordered_output_columns]
    df[reordered_output_columns].to_csv(args.output_tsv, sep='\t', header=True, index=False)
    print(f"Wrote {len(output_rows)} rows to {args.output_tsv}")


if __name__ == "__main__":
    main()
