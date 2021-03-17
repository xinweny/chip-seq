#### Packages ####
library(glue)
library(dplyr)
library(ggplot2)
library(VennDiagram)

#### Config ####
setwd("~/mrc/project/chip-seq")

#### Load data ####
bed1.only <- read.table(glue("data/K562_NHS_MAPK14_rep1_only.bed"),
                        header=FALSE, sep='\t',
                        row.names=4,
                        col.names=c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue"),
                        check.names=FALSE)
bed2.only <- read.table(glue("data/K562_HS_MAPK14_rep1_only.bed"),
                        header=FALSE, sep='\t',
                        row.names=4,
                        col.names=c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue"),
                        check.names=FALSE)
bed1bed2 <- read.table(glue("data/K562_NHS_HS_MAPK14_rep1.bed"),
                       header=FALSE, sep='\t',
                       row.names=4,
                       col.names=c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue"),
                       check.names=FALSE)

#### Venn diagram ####
dev.off()
draw.pairwise.venn(area1=nrow(bed1.only),
                   area2=nrow(bed2.only),
                   cross.area=nrow(bed1bed2),
                   scaled=FALSE,
                   category=c("NHS", "HS"))


