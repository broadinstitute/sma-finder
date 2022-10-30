## SMA Finder  

For diagnosing spinal muscular atrophy (SMA) using exome or genome sequencing data.
The tool takes 1 or more alignment files (CRAM or BAM) and reports if 
the sample(s) likely have SMA or not. It does not report carrier status or SMN2 copy number. 


### Install

```
python3 -m pip install sma-finder
```


### Usage

Example command and output:
```
/$ sma_finder --verbose --genome-version 38 --reference-fasta /ref/hg38.fa  sample1.cram

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

Full usage help text:
```

/$ sma_finder --help

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



### Output Columns

The output .tsv contains the following columns:
```
filename_prefix          =  the CRAM or BAM filename prefix 
file_type                =  "cram" or "bam"
sample_id                =  sample id from the CRAM or BAM file header (parsed from the read group metadata)  
sma_status               =  possible values are:   "has SMA",  "does not have SMA",  or "not enough coverage at SMN c.840 position"
confidence_score         =  PHRED-scaled integer confidence score that the sma_status is correct. This is similar to the PL field in GATK HaplotypeCaller genotypes.
c840_reads_with_smn1_base_C     = number of reads that have a 'C' at the c.840 position in SMN1 or SMN2  
c840_total_reads    = total number of reads overlapping the c.840 position in SMN1 plus SMN2  