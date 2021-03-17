#### Packages ####
import os
import pandas as pd
import pybedtools as pbt

#### Config ####
os.chdir("/Users/pomato/mrc/project/chip-seq/")

#### Load files ####
deseq_output = pd.read_csv("./data/GSE112379_genes_control.txt", header=0, sep='\t')
gene_bed = pd.read_csv("./data/genes_Gene.bed", header=None, names=['chr', 'start', 'end', 'name', 'score', 'strand'], sep='\t')

select_genes = list(deseq_output['geneID'])

gene_bed['name'] = gene_bed['name'].str.replace('\\.[0-9]*', '', regex=True)
filt_gene_bed = gene_bed[gene_bed['name'].isin(select_genes)]

pbt.BedTool.from_dataframe(filt_gene_bed).saveas("./data/GSE112379_genes_unchanged.bed")