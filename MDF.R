#!/usr/bin/env Rscript
library(argparser)

parser <- argparser::arg_parser( 'Produces a fasta file allowing the percentage of missing data in each site specified in -p',
                                 name = 'MDF.R')

parser <- add_argument(parser, 'input',
                       type='character',
                       nargs=1,
                       help='Path to the snpTable.tsv produced by MultiVCFAnalyzer or file.fasta. Note: if file.fasta is used as input, the flag --fasta must be included for proper parsing')
parser <- add_argument(parser, 'percentage',
                       type="integer",
                       nargs=1,
                       help='Maximum percentage of missing data for a site to be included')
parser <- add_argument(parser, 'output',
                       type='character',
                       nargs=1,
                       help='Path to output directory')
parser <- add_argument(parser, '--toExclude',
                       type = 'character',
                       help = 'Path to file containing genomes to exclude one per line',
                       default = 'None')
parser <- add_argument(parser, '--fasta',
                       type = 'logical',
                       help = 'Use when input is a fasta file',
                       flag = T, short = '-f')
argv <- parse_args(parser)

library(tidyverse)
library(seqinr)

excludingGenomes <- function(df, genometoExclude) {
  print(genometoExclude)
  toExclude <- read.delim(genometoExclude, header = F)
  df_sampleRows <- df %>%
    filter(!Genome %in% toExclude$V1)
  print(paste("Provided genomes", toExclude$V1, "have been excluded", sep = " "))
  return(df_sampleRows)
}

partialDeletionFasta <- function(input, inputdf, percentage, outputpath, genomeToExclude) {
  if( genomeToExclude == "None"){
    df <- inputdf
    print("No genomes excluded")
  } else {
    df <- excludingGenomes(inputdf, genomeToExclude)
  }
  partialDelInfo <- df %>%
    group_by(Position) %>%
    summarise(missing=sum(Call== "N"), total=n()) %>%
    mutate(percent = (missing/total)*100) %>%
    filter(percent <= percentage)
  print("Missing information calculated")
  if( input == "fasta") {
    snpTable_pd <- df %>%
      filter(Position %in% partialDelInfo$Position) %>%
      mutate(Call = as.character(Call)) %>%
      mutate(Final_Call = if_else(Call == '?'| Call == '-' , 'N' , Call)) %>%
      select(-Call) %>%
      pivot_wider(names_from = Position, values_from = Final_Call)
  } else {
    snpTable_pd <- df %>%
      filter(Position %in% partialDelInfo$Position) %>%
      mutate(Call = as.character(Call), Ref = as.character(Ref)) %>%
      mutate(Final_Call = if_else(Call == '.', Ref, Call)) %>%
      select(-Ref, -Call) %>%
      pivot_wider(names_from = Position, values_from = Final_Call)
  }


  outputTable <- paste(outputpath, percentage, "pd.tsv", sep = "_")

  write_tsv(snpTable_pd, outputTable, col_names = F)
  print("snpTable_pd saved")
  outputFasta <- paste(outputpath, percentage, "pd.fasta", sep = "_")
  system(paste("awk '{print \">\"$1; $1=\"\"; print $0}' OFS=", outputTable, ">", outputFasta))
  write_tsv(snpTable_pd, outputTable, col_names = T)
  print("Fasta_dp saved")
  }

fastatoTibble <- function(fasta) {
  fastaList <- seqinr::read.fasta(fasta, whole.header = T, forceDNAtolower = F)
  fastaM <- as_tibble(matrix(unlist(fastaList),
                             nrow = length(fastaList),
                             byrow = T),
                      .name_repair = "universal")
  fastaM$Genome <- names(fastaList)
  df_sampleRows <- fastaM %>%
    pivot_longer(-contains("Genome"), names_to = "Position",values_to = "Call")
  return(df_sampleRows)
  print("fastaToTibble Comple")
}

snpTabletoTibble <- function(snpTable) {
  df <- read.delim(snpTable, stringsAsFactors = F, check.names = F)
  df_sampleRows <- df %>%
    pivot_longer(-contains(c("Position","Ref")), names_to = "Genome", values_to = "Call")
  return(df_sampleRows)
  print("snpTableToTibble complete")
  }


inFasta=argv$fasta
print(inFasta)
print(argv$toExclude)

if( argv$fasta == T ){
  print("Fasta as input, start run with fasta mode")
  fastaT <- fastatoTibble(fasta = argv$input)
  partialDeletionFasta(input = "fasta", inputdf = fastaT, percentage = argv$percentage, outputpath = argv$output, genomeToExclude = argv$toExclude)
} else {
  print("SnpTable as input, start run with table mode")
  snpTableT <- snpTabletoTibble(snpTable = argv$input)
  partialDeletionFasta(input = "snpTable", inputdf = snpTableT, percentage = argv$percentage, outputpath = argv$output, genomeToExclude = argv$toExclude)
}
