#!/usr/env python3

"""This script computes read counts at key positions in the SMN1 & SMN2 genes and uses these to determine whether
a sample has spinal muscular atrophy (SMA).
"""

import argparse
import math
import os
import pysam
import pandas as pd
import pprint
import re
from scipy.stats import binom

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

"""The IlluminaCopyNumberCaller paper [Chen 2020] Fig 3C. considers the 2,504 unaffected individuals from the 
1kGP projects and shows that the most extreme observed ratio of SMN1 vs. SMN2 copy number is 1 to 4. Specifically, 
~10 out of 2,504 individuals (0.4%) have 4 copies of SMN2 while having only 1 copy of SMN1.  
For these individuals, we'd expect 20% of the reads that overlap the c.840 position in SMN1 and SMN2 to have the 
'C' base found in the SMN1 paralog:

      c840_reads_with_smn1_base_C / c840_total_reads ~= 0.2
"""

MAX_TOTAL_SMN_COPIES = 5


"""To differentiate individuals who have more than 0 copies of SMN1 (and so are unaffected or carriers) from 
individuals who have 0 copies of SMN1 (and so should be called as affected with SMA), we need at least 14 reads coverage 
to be certain (p < 0.05). 
"""
MIN_COVERAGE_NEEDED_TO_CALL_SMA_STATUS = 14

"""Allow for a base sequencing error rate of Q30 on the Phred scale"""
BASE_ERROR_RATE = 0.001


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


def is_zero_copies_of_smn1_more_likely_than_one_or_more_copies(n_reads_supporting_smn1, total_reads, base_error_rate):
    """Compute the likelihood of 0 copies of SMN1 vs the likelihood of 1 or more copies given the read data.

    Args:
        n_reads_supporting_smn1 (int): number of reads that support the presence of a functional SMN1 paralog
        total_reads (int): coverage estimate to use for comparison with n_reads_supporting_smn1
        base_error_rate (float): probability of a sequencing error at any given base

    Returns:
        bool: returns True if n_reads_supporting_smn1 is too large (relative to total_reads) to be treated as just a
            sequencing error. Otherwise, returns False.
    """

    n_smn1_copies_with_max_likelihood = None
    max_likelihood = 0
    for n_smn1_copies in range(0, MAX_TOTAL_SMN_COPIES + 1):
        if n_smn1_copies == 0:
            p_smn1_read = base_error_rate
            n_total = n_reads_supporting_smn1
        else:
            p_smn1_read = n_smn1_copies/MAX_TOTAL_SMN_COPIES
            n_total = total_reads

        current_likelihood = binom.pmf(n_reads_supporting_smn1, n_total, p_smn1_read)
        if current_likelihood > max_likelihood:
            max_likelihood = current_likelihood
            n_smn1_copies_with_max_likelihood = n_smn1_copies

    return n_smn1_copies_with_max_likelihood


def call_sma_status(genome_version, output_row):
    """Determines SMA status.

    Args:
        genome_version (str): sample genome version (eg. "38")
        output_row (dict): fields that will be written to the output .tsv
    """

    exon_positions = SMN_OTHER_EXON_POSITIONS_1BASED[genome_version]
    average_coverage_of_other_exons = sum(output_row[f"reads_at_{exon_name}"] for exon_name in exon_positions)
    average_coverage_of_other_exons /= len(exon_positions)

    if output_row['c840_total_reads'] >= MIN_COVERAGE_NEEDED_TO_CALL_SMA_STATUS:
        if is_zero_copies_of_smn1_more_likely_than_one_or_more_copies(
                output_row['c840_reads_with_smn1_base_C'], output_row['c840_total_reads'], BASE_ERROR_RATE):
            sma_status = "has SMA"
            sma_status_details = "has 0 copies of SMN1"
        else:
            sma_status = "doesn't have SMA"
            sma_status_details = "has 1 or more copies of SMN1"
    else:
        if not is_zero_copies_of_smn1_more_likely_than_one_or_more_copies(
                output_row['c840_reads_with_smn1_base_C'], average_coverage_of_other_exons, BASE_ERROR_RATE):
            sma_status = "may have SMA"
            sma_status_details = (
                    "has 0 copies of SMN exon 7 but has coverage of other exons of SMN, suggesting a deletion of "
                    "exon 7 in all copies of SMN1 and SMN2"
            )

        else:
            sma_status = "not enough coverage of SMN"
            sma_status_details = (
                f"average coverage of all SMN exons is < {MIN_COVERAGE_NEEDED_TO_CALL_SMA_STATUS} reads, suggesting "
                "either technical issues, or the complete deletion of both SMN1 and SMN2"
            )

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
