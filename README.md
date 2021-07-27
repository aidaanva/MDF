# MDF.R
MDF.R filters the snpTable.tsv from MultiVCFAnalyzer or a fasta alignment by removing sites that contain more than a specified percent of missing data.

## Usage
```bash
Rscript MDF.R [--] [--help] [-fasta] [--opts OPTS] input percentage output
```

## Author
Aida Andrades Valtueña (aida.andrades[at]gmail.com

## Description
```positional arguments:
  input            Path to the snpTable.tsv produced by
                   MultiVCFAnalyzer or file.fasta. Note: if file.fasta
                   is used as input, the flag --fasta must be included
                   for proper parsing
  percentage       Maximum percentage of missing data for a site to be
                   included
  output           Path to output directory

flags:
  -h, --help       show this help message and exit
  -f, --fasta      Use when input is a fasta file

optional arguments:
  -x, --opts       RDS file containing argument values
  -t, --toExclude  Path to file containing genomes to exclude one per
                   line [default: None]

```
