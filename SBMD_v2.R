#!/usr/bin/env Rscript
library(argparser)

parser <- argparser::arg_parser( 'Produces a fasta file allowing the percentage of missing data in each site specified in -p',
                                 name = 'SBMD.R')

parser <- add_argument(parser, 'input', 
                       type='character', 
                       nargs=1, 
                       help='Path to the snoTable.tsv produced by MultiVCFAnalyzer')
parser <- add_argument(parser, 'percentage', 
                       type="integer", 
                       nargs=1, 
                       help='Maximum percentage of missing data for a site to be included')
parser <- add_argument(parser, 'output', 
                       type='character', 
                       nargs=1, 
                       help='Path to output directory')
parser <- add_argument(parser, 'toExclude',
                       type = 'character',
                       nargs=1, 
                       help='File containing genomes to exclude one per line')
parser <- add_argument(parser, '-fasta',
                       type = 'logical',
                       help = 'Use when inpus is a fasta file',
                       flag = T, short = '-f')
argv <- parse_args(parser)

library(tidyverse)
library(seqinr)

#partialDeletionFasta <- function(path, percentage, outputpath, genomeToExclude) {
  df <- read.delim(path)
  toExclude <- read.delim(genomeToExclude, header = F)
  
  df_sampleRows <- df %>%
    pivot_longer(-contains(c("Position","Ref")), names_to = "Genome", values_to = "Call") %>% 
    filter(!Genome %in% toExclude$V1)
  
  partialDelInfo <- df_sampleRows %>%
    group_by(Position) %>%
    summarise(missing=sum(Call== "N"), total=n()) %>%
    mutate(percent = (missing/total)*100) %>%
    filter(percent <= percentage)
  
  snpTable_pd <- df_sampleRows %>% 
    filter(Position %in% partialDelInfo$Position) %>% 
    mutate(Call = as.character(Call), Ref = as.character(Ref)) %>% 
    mutate(Final_Call = if_else(Call == '.', Ref, Call)) %>% 
    select(-Ref, -Call) %>% 
    pivot_wider(names_from = Position, values_from = Final_Call)
  
  pd <- 100-percentage
  
  outputTable <- paste(outputpath, pd, "pd.tsv", sep = "_")
  
  write_tsv(snpTable_pd, outputTable, col_names = F)
  
  outputFasta <- paste(outputpath, pd, "pd.fasta", sep = "_")
  system(paste("awk '{print \">\"$1; $1=\"\"; print $0}' OFS=", outputTable, ">", outputFasta))
}
#fastatoTibble <- function(fasta) {
  fastaList <- read.fasta(fasta)
  fastaM <- as_tibble(matrix(unlist(fasta), nrow = length(fasta), byrow = T))
}

#snpTabletoTibble <- (snpTable) {
  df <- read.delim(path)
  
}
  
#partialDeletionFasta(argv$input, argv$percentage, argv$output, argv$toExclude)
