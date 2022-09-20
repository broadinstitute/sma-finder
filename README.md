sma_finder.py: A tool for diagnosing spinal muscular atrophy (SMA) based on exome or genome sequencing data. 

### Install

```
python3 -m pip install sma-finder
```


### Usage

```
usage: sma_finder [-h] -R REFERENCE_FASTA -g {37,38,T2T} [-o OUTPUT_TSV] [-v]
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
                        Output tsv file path
  -v, --verbose         Whether to print extra details during the run

```