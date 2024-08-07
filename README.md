## SMA Finder  

SMA Finder is a tool for diagnosing spinal muscular atrophy (SMA) from short read exome, genome, or panel sequencing data. It takes a reference sequence (FASTA) and 1 or more alignment files (CRAM or BAM) as input, evaluates reads at the c.840 position of *SMN1* and *SMN2* to detect the most common molecular causes of SMA, and then reports whether the data indicates a complete loss of functional *SMN1*, and therefore a diagnosis of SMA. 

SMA Finder has been tested on over 30,000 exome, genome, and panel sequencing samples from CMG/GREGoR rare disease cohorts as well as 200,000 exomes and genomes from the UK Biobank (UKBB). In these tests, SMA Finder's false positive rate was less than 1 in 200,000, its true positive rate was 28/28 (100%), and its positive predictive value (PPV) was 97%. 

*Limitations:*  
- does not report SMA carrier status or *SMN1/SMN2* copy numbers  
- does not detect the ~5% of cases caused by *SMN1* loss-of-function mutations that do not involve the c.840 position  
- requires at least 14 reads to overlap the c.840 position in *SMN1* plus *SMN2* in order to make a call  
- was developed and tested on Illumina short read sequencing data generated using DNA extrated from whole blood and aligned using the BWA aligner. Performance on data from other sequencing technologies, sample types, and alignment pipelines is unknown. 

### Publication

For more information about SMA Finder, see:

*Weisburd B, Sharma R, Pata V, et al. Detecting missed diagnoses of spinal muscular atrophy in genome, exome, and panel sequencing datasets. Preprint. medRxiv. 2024;2024.02.11.24302646. Published 2024 Feb 27. [doi:10.1101/2024.02.11.24302646](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10889006/)*

### Install

To install the latest version of SMA Finder, run:
```
python3 -m pip install -U sma-finder
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
---
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
<br />  

---

  
### Combining results from multiple samples

After running SMA Finder on many samples, it's often useful to combine the per-sample output tables into
a single table. One way to do this is with the following shell command:

```
cd <directory with multiple SMA Finder output tsvs>
combined_table_filename=combined_results.tsv
head -n 1 $(ls *.sma_finder_results.tsv | head -n 1) > ${combined_table_filename}   # get table header from the 1st table 
for i in *.sma_finder_results.tsv; do
    tail -n +2 $i >> ${combined_table_filename}    # concatenate all tables
done
```

### Plotting combined results

A scatter plot summarizing read counts from many samples can be generated using the `plot_SMN1_SMN2_scatter` command:

```
python3 plot_SMN1_SMN2_scatter.py --format svg --format png ${combined_table_filename}
```

It generates plots like this one which is based on 16,626 exomes that include neuromuscular disease cohorts:

<img width="799" alt="image" src="https://github.com/broadinstitute/sma-finder/assets/6240170/d097a231-9b66-445b-b53c-84abdb9887d0">

---
### Poster from SVAR22

This poster on SMA Finder was presented at the [SVAR22](https://www.grahamerwin.org/svar-conference) conference:

<img src="https://github.com/broadinstitute/sma_finder/raw/main/docs/SMA_poster_SVAR22.png" />

