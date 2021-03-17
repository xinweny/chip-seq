#### Packages ####
library(tidyverse)

setwd("/Users/Pomato/mrc/project/chip-seq")

#### Functions ####
read_deeptools_table <- function(file) {
  
  n <- max(count.fields(file, sep = '\t'), na.rm = TRUE)
  x <- readLines(file)
  
  .splitvar <- function(x, sep, n) {
    var <- unlist(strsplit(x, split = sep))
    length(var) <- n
    return(var)
  }
  
  x <- do.call(cbind, lapply(x, .splitvar, sep = '\t', n = n))
  x <- apply(x, 1, paste, collapse = '\t')
  plot_table <- na.omit(read.csv(text = x, sep = '\t')[-1,])  # Remove first row with "gene" label
  
  return(plot_table)
}

#### Load data ####
table <- read_deeptools_table('./data/GSE112379_K562_RPB1_HS_vs_NHS_profilePlot_scale-regions.tab')

