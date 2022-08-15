#!/usr/env python3

"""This script counts the number of reads that overlap exons of SMN1 and SMN2 in the input CRAM file."""

import argparse
import os
import pysam
import pandas as pd
import re

SMN_COORDINATES = {
    "SMN_CHROMOSOME": {
        "37": "5",
        "38": "chr5",
        "T2T": "5",
    },
    "SMN1_C840_POSITION_1BASED": {
        "37": 70247773,
        "38": 70951946,
        "T2T": 71408734,
    },
    "SMN2_C840_POSITION_1BASED": {
        "37": 69372353,
        "38": 70076526,
        "T2T": 70810812,
    }
}

VALID_GENOME_VERSIONS = set(SMN_COORDINATES["SMN_CHROMOSOME"].keys())


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
    """Count the number of A, C, G, or T found at the given genomic position within the given read data.

    Args:
        alignment_file (pysam.AlignmentFile): A pysam.AlignmentFile object representing a BAM or CRAM file.
        chrom (str): Chromosome of the genomic position.
        pos_1based (int): 1-based genomic position where to count nucleotides.

    Return:
        dict: The keys are nucleotides "A", "C", "G", "T", and the values are counters representing the number of times
            a read contained that nucleotide at the given position.
    """

    pos_0based = pos_1based - 1
    nucleotide_counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    for pileup_column in alignment_file.pileup(region=f"{chrom}:{pos_0based}-{pos_1based}", truncate=True):
        if pileup_column.pos < pos_0based:
            continue

        if pileup_column.pos != pos_0based:
            raise ValueError(f"Unexpected pileup position: {chrom}:{pileup_column.pos}. Expecting {chrom}:{pos_0based}")

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


def main():
    args = parse_args()

    if args.verbose:
        print(f"Processing {len(args.cram_or_bam_path)} CRAM files")

    output_records = []
    for cram_or_bam_path in args.cram_or_bam_path:
        print(f"==> {cram_or_bam_path}")
        alignment_filename = os.path.basename(cram_or_bam_path)
        sample_id = re.sub("(.cram|.bam)$", "", alignment_filename)

        chrom = SMN_COORDINATES["SMN_CHROMOSOME"][args.genome_version]
        smn1_c840_position = SMN_COORDINATES["SMN1_C840_POSITION_1BASED"][args.genome_version]
        smn2_c840_position = SMN_COORDINATES["SMN2_C840_POSITION_1BASED"][args.genome_version]

        with pysam.AlignmentFile(cram_or_bam_path, 'rc', reference_filename=args.reference_fasta) as alignment_file:
            smn2_c840_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, smn2_c840_position)
            smn1_c840_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, smn1_c840_position)

        smn1_c840_total_reads = sum(smn1_c840_nucleotide_counts.values())
        smn2_c840_total_reads = sum(smn2_c840_nucleotide_counts.values())
        smn_c840_total_reads = smn1_c840_total_reads + smn2_c840_total_reads
        if args.verbose:
            print(f"{smn_c840_total_reads:15,d} total reads in {cram_or_bam_path}")

        output_records.append({
            "sample_id": sample_id,
            "smn_c840_reads_with_C": smn1_c840_nucleotide_counts["C"] + smn2_c840_nucleotide_counts["C"],
            "smn_c840_total_reads": smn_c840_total_reads,
            "smn1_c840_reads_with_C": smn1_c840_nucleotide_counts["C"],
            "smn2_c840_reads_with_C": smn2_c840_nucleotide_counts["C"],
        })

    df = pd.DataFrame(output_records)
    df.to_csv(args.output_tsv, sep='\t', header=True, index=False)
    print(f"Wrote {len(output_records)} rows to {args.output_tsv}")


if __name__ == "__main__":
    main()
