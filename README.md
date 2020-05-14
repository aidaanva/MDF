# MDF.R
MDF.R filters a fasta alignment by removing sites that contain more than a specified percent of missing data.

## Usage
```bash
Rscript MDF.R [--] [--help] [-fasta] [--opts OPTS] input percentage output toExclude
```

## Author
Aida Andrades Valtueña (aida.andrades[at]gmail.com

## Description
```bash
positional arguments:
  input       Path to the snpTable.tsv produced by MultiVCFAnalyzer
  percentage  Maximum percentage of missing data for a site to be
              included
  output      Path to output directory
  toExclude   File containing genomes to exclude one per line

flags:
  -h, --help  show this help message and exit
  -f, -fasta  Use when inpus is a fasta file

optional arguments:
  -x, --opts  RDS file containing argument values

```
