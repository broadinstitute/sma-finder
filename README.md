## SMA Finder  

A tool for diagnosing spinal muscular atrophy (SMA) using exome or genome sequencing data.
It takes 1 or more alignment files (CRAM or BAM) as input and reports whether 
the sample(s) have complete loss of functional *SMN1*.

Testing on 13 positive controls and 10,434 negative controls showed 100% sensitivity and specificity (details below). 

*Limitations:* SMA Finder doesn't report SMA carrier status or *SMN2* copy number. Also, it does not detect the ~5% of cases caused by unusual *SMN1* loss-of-function mutations that do not involve the c.840 position. Finally, although SMA Finder may also work for RNA-seq, targeted sequencing, or long-read alignment files, we currently lack the positive controls needed to evaluate this. If you have such samples, or would like to collaborate in other ways, please contact weisburd at broadinstitute dot org.


### Install

```
python3 -m pip install sma-finder
```

### Example

Example command:
```
sma_finder --verbose --genome-version 38 --reference-fasta /ref/hg38.fa  sample1.cram
```
Command output:
```
Input args:
    --reference-fasta: /ref/hg38.fa
    --genome-version: 38
    --output-tsv: sample1.sma_finder_results.tsv
    CRAMS or BAMS: sample1.cram
----
Output row #1:
        filename_prefix                     sample1
        file_type                           cram
        sample_id                           s1
        sma_status                          has SMA
        confidence_score                    168
        c840_reads_with_smn1_base_C         0
        c840_total_reads                    174
Wrote 1 rows to sample1.sma_finder_results.tsv        
```

### Usage

Usage help text:

```
sma_finder --help

usage: sma_finder [-h] -R REFERENCE_FASTA -g {37,38,T2T} [-o OUTPUT_TSV]
                     [-v]
                     cram_or_bam_path [cram_or_bam_path ...]

positional arguments:
  cram_or_bam_path      One or more CRAM or BAM file paths

optional arguments:
  -h, --help            show this help message and exit
  -R REFERENCE_FASTA, --reference-fasta REFERENCE_FASTA
                        Reference genome FASTA file path
  -g {37,38,T2T}, --genome-version {37,38,T2T}
                        Reference genome version
  -o OUTPUT_TSV, --output-tsv OUTPUT_TSV
                        Optional output tsv file path
  -v, --verbose         Whether to print extra details during the run
```

---

### Output

The output .tsv contains one row per input CRAM or BAM file and has the following columns:

<table>
    <tr>
        <td><b>filename_prefix</b></td>
        <td>CRAM or BAM filename prefix. If the input file is <i>/path/sample1.cram</i> this would be <i>"sample1"</i>.</td>
    </tr>
    <tr>
        <td><b>file_type</b></td>
        <td><i>"cram"</i> or <i>"bam"</i></td>
    </tr>
    <tr>
        <td><b>sample_id</b></td>
        <td>sample id from the CRAM or BAM file header (parsed from the read group metadata)</td>
    </tr>
    <tr>
        <td><b>sma_status</b></td>
        <td>possible values are:<br> 
            <i>"has SMA"</i><br>
            <i>"does not have SMA"</i><br>
            <i>"not enough coverage at SMN c.840 position"</i><br>
        </td>
    <tr>
        <td><b>confidence_score</b></td>
        <td>PHRED-scaled integer score measuring the level of confidence that the sma_status is correct. The bigger the score, the higher the confidence. It is calculated in a similar way to the PL field in GATK HaplotypeCaller genotypes.</td>
    <tr>
        <td><b>c840_reads_with_smn1_base_C</b></td>
        <td>number of reads that have a 'C' nucleotide at the c.840 position in SMN1 plus SMN2</td> 
    <tr>
        <td><b>c840_total_reads</b></td>
        <td>total number of reads overlapping the c.840 position in SMN1 plus SMN2</td>  
    </tr>
</table>

---
### Details

This poster on SMA Finder was presented at the [SVAR22](https://www.grahamerwin.org/svar-conference) conference:

<img src="https://github.com/broadinstitute/sma_finder/raw/main/docs/SMA_poster_SVAR22.png" />
