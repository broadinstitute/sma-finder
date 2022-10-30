#!/usr/env python3

"""Spinal muscular atrophy (SMA) is a recessive genetic disease caused by the SMN gene. SMN consists of 2
nearly-identical paralogs: SMN1 and SMN2. Both are located on chromosome 5 and differ from each other at only 2 exonic
positions: c.840 and c.1124. SMN1 contains a 'C' nucleotide at the c.840 position while SMN2 contains a 'T'.
This key difference causes skipping of exon 7 in most transcripts of SMN2 and explains why SMN1 is the main source of
functional SMN protein. An individual's genome can have 0 or more copies of SMN1 and 0 or more copies of SMN2.
If an individual has exactly 0 copies of functional SMN1, then they develop SMA.
This script determines whether a sample is positive for SMA by examining read counts at the key c.840 position in
SMN1 & SMN2.
"""

import argparse
import hashlib
import math
import os
import pysam
import pandas as pd
import re
from scipy.stats import binom

"""Output column names that will be populated for each sample in the output table"""
OUTPUT_COLUMNS = [
    "filename_prefix",
    "file_type",
    "sample_id",
    "sma_status",
    "confidence_score",
    "c840_reads_with_smn1_base_C",
    "c840_total_reads",
]


"""Each key is a reference genome version, and each value is a 5-tuple with the genomic coordinates and reference base 
at the c.840 positions in SMN1 and SMN2: 

(
    chromosome, 
    1-based coordinates of the c.840 position in SMN1, 
    SMN1 reference base, 
    1-based coordinates of the c.840 position in SMN2, 
    SMN2 reference base
)  
"""
SMN_C840_POSITION_1BASED = {
    "37": ("5", 70247773, "C", 69372353, "T"),
    "38": ("chr5", 70951946, "C", 70076526, "T"),
    "T2T": ("5", 71408734, "C",  70810812, "T"),
}

"""MAX_TOTAL_SMN_COPIES represents the largest number of total SMN copies we'd expect to see when an individual has only 
1 copy of SMN1.    

The IlluminaCopyNumberCaller paper [Chen 2020] Fig 3C. considers 2,504 unaffected individuals from the 1kGP project and 
shows that the most extreme observed ratio of SMN1 vs. SMN2 copy number is 1 to 4. Specifically, ~10 out of 2,504 
individuals (0.4%) have 4 copies of SMN2 while having only 1 copy of SMN1.  
"""
MAX_TOTAL_SMN_COPIES = 5


"""To confidently distinguish individuals who have who have 0 copies of SMN1 (and so are affected with SMA) from 
individuals with 1 or more copies of SMN1 (who are unaffected or carriers), we want there to be at least 14 reads 
overlapping the c.840 position. This is because, with a null hypothesis of 1 copy of SMN1 and 4 copies of SMN2, the 
binomial(r=0, n=14, p=1/5) yields p < 0.05.
"""
MIN_COVERAGE_NEEDED_TO_CALL_SMA_STATUS = 14

"""Assume a base sequencing error of Q30 (Phred scale) which = 0.001. Divide by 3 since we only care about errors 
that convert the reference base to 1 of the 3 possible other bases (ie. convert the reference 'T' at the c.840 position 
in SMN2 to a 'C' which matches the c.840 position in SMN1). This is similar to equation 6 in the [Poplin 2018] 
HaplotypeCaller paper - in the section that describes the reference model. 
"""
BASE_ERROR_RATE = 0.001 / 3

"""A minimum value for log-likelihood. This is used when very small likelihood values (eg. 10**-500) are rounded down
to zero and lead to domain errors in the log calculation.   
"""
MINIMUM_LOG_LIKELIHOOD = -300


def parse_args():
    """Define and then parse command-line args"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference genome FASTA file path")
    parser.add_argument("-g", "--genome-version", required=True, choices=sorted(SMN_C840_POSITION_1BASED.keys()),
                        help="Reference genome version")
    parser.add_argument("-o", "--output-tsv", help="Optional output tsv file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Whether to print extra details during the run")
    parser.add_argument("cram_or_bam_path", nargs="+", help="One or more CRAM or BAM file paths")

    args = parser.parse_args()

    # make sure the input files exist
    for path in [args.reference_fasta] + args.cram_or_bam_path:
        if not os.path.isfile(path):
            parser.error(f"File not found: {path}")

    # define the output_tsv if it wasn't specified
    if not args.output_tsv:
        if len(args.cram_or_bam_path) > 1:
            input_paths_hash = hashlib.sha256("|".join(args.cram_or_bam_path).encode('utf-8')).hexdigest().upper()
            args.output_tsv = f"sma_finder_results.{len(args.cram_or_bam_path)}_samples.{input_paths_hash[:10]}.tsv"
        else:
            filename_prefix, _ = get_filename_prefix_and_file_type(args.cram_or_bam_path[0])
            args.output_tsv = f"{filename_prefix}.sma_finder_results.tsv"

    if args.verbose:
        print("Input args:")
        print(f"    --reference-fasta: {os.path.abspath(args.reference_fasta)}")
        print(f"    --genome-version: {args.genome_version}")
        print(f"    --output-tsv: {os.path.abspath(args.output_tsv)}")
        print(f"    CRAMS or BAMS:", ", ".join(map(os.path.abspath, args.cram_or_bam_path)))

    return args


def count_nucleotides_at_position(alignment_file, chrom, pos_1based):
    """Count the number of A, C, G, and T's found at the given genomic position within the given read data.

    Args:
        alignment_file (pysam.AlignmentFile): The pysam AlignmentFile object representing the input BAM or CRAM file.
        chrom (str): Chromosome of the genomic position.
        pos_1based (int): 1-based genomic position where to count nucleotides.

    Return:
        dict: The keys are nucleotides "A", "C", "G", "T", and the values are counts representing the number of times
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

        # iterate over the reads in the pileup
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


def get_filename_prefix_and_file_type(cram_or_bam_path):
    """Returns the filename prefix and file type suffix.

    Args:
        cram_or_bam_path (str): Input CRAM or BAM path.
        output_row (dict): dictionary of fields to be written to the .tsv output row for the current sample
    """
    filename_pattern_match = re.match("(.*)[.](cram|bam)$", os.path.basename(cram_or_bam_path))
    if not filename_pattern_match:
        raise ValueError(f"File path doesn't end with have '.cram' or '.bam': {cram_or_bam_path}")

    filename_prefix = filename_pattern_match.group(1)
    file_type = filename_pattern_match.group(2)

    return filename_prefix, file_type


def set_sample_id(alignment_file, output_row):
    """Try reading the sample_id from the alignment file header by looking for a read group (@RG) with a sample (SM)
    field.

    Args:
        alignment_file (pysam.AlignmentFile): The pysam AlignmentFile object representing the input BAM or CRAM file.
        output_row (dict): dictionary of fields to be written to the .tsv output row for the current sample
    """
    output_row["sample_id"] = ""
    for read_group in alignment_file.header.get("RG", []):
        if "SM" in read_group:
            output_row["sample_id"] = read_group["SM"]
            break


def is_zero_copies_of_smn1_more_likely_than_one_or_more_copies(n_reads_supporting_smn1, total_reads, base_error_rate):
    """Compute the likelihood of the read data given 0 copies of SMN1 vs the likelihood given 1 or more copies.

    Args:
        n_reads_supporting_smn1 (int): Number of reads that support the presence of a functional SMN1 paralog
        total_reads (int): Coverage estimate to use for comparison with n_reads_supporting_smn1
        base_error_rate (float): Probability of a sequencing error at any given base

    Returns:
        bool: True if the likelihood of the given read data is greater if you assume 0 copies of SMN1 than if you
            assume 1 or more copies of SMN1.
        float: A PHRED-scale confidence score
    """
    if total_reads == 0:
        raise ValueError("total_reads == 0")

    likelihood_given_zero_copies = binom.pmf(n_reads_supporting_smn1, total_reads, base_error_rate)

    likelihood_given_one_or_more_copies = 0
    for n_smn1_copies in range(1, MAX_TOTAL_SMN_COPIES + 1):
        p_smn1_read = n_smn1_copies/MAX_TOTAL_SMN_COPIES
        likelihood_given_one_or_more_copies = max(
            likelihood_given_one_or_more_copies,
            binom.pmf(n_reads_supporting_smn1, total_reads, p_smn1_read))

    if likelihood_given_zero_copies == 0 and likelihood_given_one_or_more_copies == 0:
        raise ValueError(f"likelihood_given_zero_copies == 0 and likelihood_given_one_or_more_copies == 0 when "
                         f"n_reads_supporting_smn1={n_reads_supporting_smn1}, total_reads={total_reads}, and "
                         f"base_error_rate={base_error_rate}")

    log_likelihood_given_zero_copies = MINIMUM_LOG_LIKELIHOOD
    if likelihood_given_zero_copies > 0:
        log_likelihood_given_zero_copies = max(math.log10(likelihood_given_zero_copies), MINIMUM_LOG_LIKELIHOOD)

    log_likelihood_given_one_or_more_copies = MINIMUM_LOG_LIKELIHOOD
    if likelihood_given_one_or_more_copies > 0:
        log_likelihood_given_one_or_more_copies = max(math.log10(likelihood_given_one_or_more_copies), MINIMUM_LOG_LIKELIHOOD)

    phred_scaled_confidence_score = int(10 * abs(
        log_likelihood_given_zero_copies - log_likelihood_given_one_or_more_copies
    ))
    return likelihood_given_zero_copies > likelihood_given_one_or_more_copies, phred_scaled_confidence_score


def call_sma_status(output_row):
    """Compute SMA status and set 'sma_status' and 'confidence_score' fields in the output_row.

    Args:
        output_row (dict): dictionary of fields to be written to the .tsv output row for the current sample
    """

    if output_row['c840_total_reads'] >= MIN_COVERAGE_NEEDED_TO_CALL_SMA_STATUS:
        zero_copies_more_likely, confidence_score = is_zero_copies_of_smn1_more_likely_than_one_or_more_copies(
            output_row['c840_reads_with_smn1_base_C'], output_row['c840_total_reads'], BASE_ERROR_RATE)
        if zero_copies_more_likely:
            sma_status = "has SMA"
        else:
            sma_status = "does not have SMA"
    else:
        sma_status = "not enough coverage at SMN c.840 position"
        confidence_score = 0

    output_row["sma_status"] = sma_status
    output_row["confidence_score"] = confidence_score


def main():
    args = parse_args()

    chrom, c840_position_in_smn1, smn1_base, c840_position_in_smn2, _ = SMN_C840_POSITION_1BASED[args.genome_version]

    # process the input BAM or CRAM files
    output_rows = []
    for cram_or_bam_path in args.cram_or_bam_path:
        output_row = {}
        output_row["filename_prefix"], output_row["file_type"] = get_filename_prefix_and_file_type(cram_or_bam_path)

        with pysam.AlignmentFile(cram_or_bam_path, 'rc', reference_filename=args.reference_fasta) as alignment_file:
            set_sample_id(alignment_file, output_row)
            smn1_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, c840_position_in_smn1)
            smn2_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, c840_position_in_smn2)

        output_row.update({
            f"c840_reads_with_smn1_base_{smn1_base}": smn1_nucleotide_counts[smn1_base] + smn2_nucleotide_counts[smn1_base],
            f"c840_total_reads": sum(smn1_nucleotide_counts.values()) + sum(smn2_nucleotide_counts.values()),
        })

        call_sma_status(output_row)

        output_rows.append(output_row)

    # write results to .tsv
    df = pd.DataFrame(output_rows)

    if args.verbose:
        for i, (_, row) in enumerate(df.iterrows()):
            print("----")
            print(f"Output row #{i+1}:")
            for column in OUTPUT_COLUMNS:
                print(f"        {column:35s} {row[column]}")

    df[OUTPUT_COLUMNS].to_csv(args.output_tsv, sep='\t', header=True, index=False)
    print(f"Wrote {len(output_rows)} rows to {os.path.abspath(args.output_tsv)}")


if __name__ == "__main__":
    main()
