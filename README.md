## SMA Finder  

SMA Finder is a tool for diagnosing spinal muscular atrophy (SMA) using exome or genome sequencing data.  
It takes a reference sequence (FASTA) and 1 or more alignment files (CRAM or BAM) as input, evaluates reads that overlap the 
c.840 position in *SMN1* and *SMN2*, and then reports whether the input sample(s) have complete loss of functional *SMN1*.

Running this tool on 13 positive controls and 10,434 negative controls showed 100% sensitivity and specificity (details below). 

*Limitations:*  
SMA Finder doesn't report SMA carrier status or *SMN2* copy number.  Also, it does not detect the ~5% of cases caused by unusual *SMN1* loss-of-function mutations that do not involve the c.840 position. 


### Install

```
python3 -m pip install sma-finder
```

### Example

Example command:
```
sma_finder --verbose --hg38-reference-fasta /ref/hg38.fa  sample1.cram
```
Command output:
```
Input args:
    --hg38-reference-fasta: /ref/hg38.fa
    --output-tsv: sample1.sma_finder_results.tsv
    CRAMS or BAMS: sample1.cram
----
Output row #1:
        filename_prefix                     sample1
        file_type                           cram
        genome_version                      hg38
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

usage: sma_finder.py [-h] [--hg37-reference-fasta HG37_REFERENCE_FASTA]
                     [--hg38-reference-fasta HG38_REFERENCE_FASTA]
                     [--t2t-reference-fasta T2T_REFERENCE_FASTA]
                     [-o OUTPUT_TSV] [-v]
                     cram_or_bam_path [cram_or_bam_path ...]

positional arguments:
  cram_or_bam_path      One or more CRAM or BAM file paths

optional arguments:
  -h, --help            show this help message and exit
  --hg37-reference-fasta HG37_REFERENCE_FASTA
                        HG37 reference genome FASTA path. This should be
                        specified if the input bam or cram is aligned to HG37.
  --hg38-reference-fasta HG38_REFERENCE_FASTA
                        HG38 reference genome FASTA path. This should be
                        specified if the input bam or cram is aligned to HG38.
  --t2t-reference-fasta T2T_REFERENCE_FASTA
                        T2T reference genome FASTA path. This should be
                        specified if the input bam or cram is aligned to the
                        CHM13 telomere-to-telomere benchmark.
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
        <td><b>genome_version</b></td>
        <td><i>"hg37"</i>, <i>"hg38"</i>, or <i>"t2t"</i></td>
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
