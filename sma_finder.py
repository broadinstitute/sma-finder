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
import re
from scipy.stats import binom
import sys

"""Output column names that will be populated for each sample in the output table"""
OUTPUT_COLUMNS = [
    "filename_prefix",
    "file_type",
    "genome_version",
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
    #"38": ("chr5_KI270897v1_alt", 500378, "C", 301867, "T"),
    #"38": ("chr5_GL339449v2_alt", 458845, "G", 458845, "T"),
    "t2t":  ("chr5", 71408734, "C",  70810812, "T"),
}

# NOTE: in GRCh38, there are alternative contigs with SMN1 and SMN2 sequences:
#   SMN1 @ chr5_KI270897v1_alt:473491-501447  with the c.840 C base @ chr5_KI270897v1_alt:500378
#   SMN2 @ chr5_KI270897v1_alt:274951-303848  with the c.840 T base @ chr5_KI270897v1_alt:301867
#   reverse complement of SMN1 @ chr5_GL339449v2_alt:456849-485731  with the c.840 G base @ chr5_GL339449v2_alt:458845
# Since, by default, bwa mem performs ALT-aware alignment, we don't expect to see primary alignments at these alternative
# contig sequences (see https://github.com/lh3/bwa/blob/master/README-alt.md#step-1-bwa-mem-mapping). Checking ~20,000
# exomes and ~5,000 genomes from cohorts @ the Broad Institute Center for Mendelian Genomics confirmed this.
# In all samples, there were 0 reads that aligned to the c.840 positions in the alternative contig sequences of SMN.

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
    parser.add_argument("--hg37-reference-fasta", help="HG37 reference genome FASTA path. This should be specified if "
                                                      "the input bam or cram is aligned to HG37.")
    parser.add_argument("--hg38-reference-fasta", help="HG38 reference genome FASTA path. This should be specified if "
                                                      "the input bam or cram is aligned to HG38.")
    parser.add_argument("--t2t-reference-fasta", help="T2T reference genome FASTA path. This should be specified if "
                                                     "the input bam or cram is aligned to the CHM13 "
                                                     "telomere-to-telomere benchmark.")
    parser.add_argument("-o", "--output-tsv", help="Optional output tsv file path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Whether to print extra details during the run")
    parser.add_argument("cram_or_bam_path", nargs="+", help="One or more CRAM or BAM file paths")

    args = parser.parse_args()

    # make sure the input files exist
    fasta_paths = {
        "37": args.hg37_reference_fasta,
        "38": args.hg38_reference_fasta,
        "t2t": args.t2t_reference_fasta,
    }
    fasta_paths = {label: path for label, path in fasta_paths.items() if path}

    if not fasta_paths:
        parser.error(f"Reference gnome fasta path not specified")

    for path in list(fasta_paths.values()) + args.cram_or_bam_path:
        if not os.path.isfile(os.path.expanduser(path)):
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
        for label, path in fasta_paths.items():
            label = f"hg{label}" if label != "t2t" else label
            print(f"    --{label}-reference-fasta: {os.path.abspath(path)}")
        print(f"    --output-tsv: {os.path.abspath(args.output_tsv)}")
        print(f"    CRAMS or BAMS:", ", ".join(map(os.path.abspath, args.cram_or_bam_path)))
        print("----")

    return args, fasta_paths


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
    args, fasta_paths = parse_args()

    # process the input BAM or CRAM files
    output_rows = []
    for cram_or_bam_path in args.cram_or_bam_path:
        for genome_version, reference_fasta_path in fasta_paths.items():
            chrom, c840_position_in_smn1, smn1_base, c840_position_in_smn2, _ = SMN_C840_POSITION_1BASED[genome_version]

            filename_prefix, file_type = get_filename_prefix_and_file_type(cram_or_bam_path)

            genome_version_label = f"hg{genome_version}" if genome_version != "t2t" else genome_version
            output_row = {
                "filename_prefix": filename_prefix,
                "file_type": file_type,
                "genome_version": genome_version_label,
            }

            try:
                with pysam.AlignmentFile(
                        cram_or_bam_path, 'rc', reference_filename=reference_fasta_path) as alignment_file:
                    set_sample_id(alignment_file, output_row)
                    smn1_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, c840_position_in_smn1)
                    smn2_nucleotide_counts = count_nucleotides_at_position(alignment_file, chrom, c840_position_in_smn2)
            except ValueError as e:
                print(f"ERROR: unable to get read counts from {cram_or_bam_path} for {genome_version_label}: {e}. "
                      f"Skipping...")
                continue
            else:
                if len(fasta_paths) > 1:
                    print(f"Retrieved read counts from {cram_or_bam_path} for {genome_version_label}")

            output_row.update({
                f"c840_reads_with_smn1_base_{smn1_base}": smn1_nucleotide_counts[smn1_base] +
                                                          smn2_nucleotide_counts[smn1_base],
                f"c840_total_reads": sum(smn1_nucleotide_counts.values()) +
                                     sum(smn2_nucleotide_counts.values()),
            })

            call_sma_status(output_row)

            output_rows.append(output_row)

    if args.verbose:
        print("----")

    # write results to .tsv
    with open(args.output_tsv, "wt") as f:
        f.write("\t".join(OUTPUT_COLUMNS) + "\n")

        for i, row in enumerate(output_rows):
            f.write("\t".join([str(row[c]) for c in OUTPUT_COLUMNS]) + "\n")

            if args.verbose:
                print(f"Output row #{i+1}:")
                for column in OUTPUT_COLUMNS:
                    print(f"        {column:35s} {row[column]}")

    print(f"Wrote {len(output_rows)} rows to {os.path.abspath(args.output_tsv)}")

    if not output_rows:
        print(f"ERROR: no output rows")
        sys.exit(1)


if __name__ == "__main__":
    main()
