#!/usr/env python3

"""This script reports read counts in SMN1 & SMN2 that can be used to check if the given WES or WGS sample has SMA."""

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


def main():
    args = parse_args()

    output_rows = []
    for cram_or_bam_path in args.cram_or_bam_path:
        filename_pattern_match = re.match("(.*)[.](cram|bam)$", os.path.basename(cram_or_bam_path))
        if not filename_pattern_match:
            raise ValueError(f"File path doesn't end with have '.cram' or '.bam': {cram_or_bam_path}")

        output_row = {
            "filename": filename_pattern_match.group(1),
            "file_type": filename_pattern_match.group(2),
        }

        alignment_file = pysam.AlignmentFile(cram_or_bam_path, 'rc', reference_filename=args.reference_fasta)

        # determine the sample_id
        output_row["sample_id"] = ""
        for read_group in alignment_file.header.get("RG", []):
            if "SM" in read_group:
                output_row["sample_id"] = read_group["SM"]
                break

        # count reads at differing bases
        chrom = SMN_CHROMOSOME[args.genome_version]
        smn_differing_positions = SMN_DIFFERING_POSITION_1BASED[args.genome_version]
        for key, smn1_base in [("c840", "C"), ("c1124", "G")]:
            differing_position = smn_differing_positions[key]
            smn1_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, differing_position["SMN1"])
            smn2_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, differing_position["SMN2"])
            output_row.update({
                f"{key}_reads_with_smn1_base_{smn1_base}": smn1_nucleotide_counts[smn1_base] + smn2_nucleotide_counts[smn1_base],
                f"{key}_total_reads": sum(smn1_nucleotide_counts.values()) + sum(smn2_nucleotide_counts.values()),
            })

        # count reads at other exons
        for exon_name, exon_positions in SMN_OTHER_EXON_POSITIONS_1BASED[args.genome_version].items():
            smn1_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, exon_positions["SMN1"])
            smn2_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, exon_positions["SMN2"])
            output_row.update({
                f"reads_at_{exon_name}": sum(smn1_nucleotide_counts.values()) + sum(smn2_nucleotide_counts.values()),
            })

        # record output row
        output_rows.append(output_row)
        if args.verbose:
            print("Output row:")
            pprint.pprint(output_row)

        alignment_file.close()

    # write results to .tsv
    df = pd.DataFrame(output_rows)
    df.to_csv(args.output_tsv, sep='\t', header=True, index=False)
    print(f"Wrote {len(output_rows)} rows to {args.output_tsv}")


if __name__ == "__main__":
    main()
