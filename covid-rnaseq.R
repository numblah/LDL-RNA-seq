#load packages
#tidy verse packages 
library(readr)
library(dplyr)
library(magrittr)
library(tximport)
library(DESeq2)
library(ggplot2)
library('plotly')
library(tibble)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
#could have just laoded library(tidyverse)

#creates a tibble an extension of a dataframe
sample_table = read_csv('SraRunTable.txt') %>% 
  select('Sample Name',source_name, treatment,
         Cell_Line, Cell_type, time_point) %>% 
  slice(seq(1,48, by=4)) # or could use unique(sample_table)
# we need to turn this dataframe into a sample table (it is now a run table)

??tximport

#gives a vector from  a tibble column

sample_files=paste0(pull(sample_table, `Sample Name`), '/quant.sf') 

#give it names (otherwise it would be just numbered)

names(sample_files) = pull (sample_table, `Sample Name`)

#names are used to give column headers so we can track the columns 
#in the resulting data 


# for tx2gene we will need to go to linux
#find the gencode file 
# ensmbl transcript ID and corresponding ensmbl gene ID
# make two files for each
#same lengths with IDs in the same order
#3 commands
# grep -P -o 'ENST\d{11}' FILENAME > enst.txt
# grep -P -o 'ENSG\d{11}' FILENAME > ensg.txt
# third:  use paste which merges lines of files
# paste -d ',' enst.txt ensg.txt > Share/covid-rnaseq-gene_map.csv
# end now we import the file

gene_map = read_csv('gene_map.csv', col_names=c('enstid', 'ensgid'))

#come back from 
count_data= tximport(files = sample_files,
         type = 'salmon',
         tx2gene = gene_map,
         ignoreTxVersion = TRUE)


#count_data= tximport(files = sample_files,
#                     type = 'salmon',
#                     tx2gene = gene_map,
#                     ignoreTxVersion = TRUE)

#gives a list of differentially expressed genes (not transcripts)
#this makes it easier


# data type - Large list
#abundance 

count_data['counts']
count_data$counts[1:6,]


# deseq wants:
# sample_table processed :
# turn it into a data frame
# remove spaces from file names
# change statement to a factor (one things or another)

sample_table = as.data.frame(sample_table)
colnames(sample_table)[1] = 'Sample'
conditions =c('mock_nhbe', 'infected_nhbe',
              'mock_a549', 'infected_a549')

conditions = rep(conditions, each=3)
conditions = factor(conditions)
sample_table[['conditions']] = conditions

# ~ a given y ~ X aka Y give X
# DESeq will thus generate a dataset looking at counts based on conditions

deseq_dataset = DESeqDataSetFromTximport(txi = count_data,
                                         colData= sample_table,
                                         design =~conditions)
# lets look at it
counts(deseq_dataset)[1:6, 1:3]

#is this enough prep work?
#no we need to normalize, otherwise we couldn't tell if diffs are real
#read map to a gene in proportion to
#1 gene length
#2gene expression ----
#3library depth (how many sequencing reads have there been for each sample)
#need to normalize for 3 (length is not bad because both things compared will#
# have the same length)

#DESeq TMM trimmed medium of means----

#per row geometric mean of counts
#per row ratio of counts to geometric mean
#median of these ratios for a sample (column) = 'scaling factor'
#raw counts/scaling factor = normalized counts

#

deseq_dataset = estimateSizeFactors(deseq_dataset)
normalizationFactors(deseq_dataset) # actually includes gene-gene normalization
counts(deseq_dataset, normalized=TRUE)[1:6, 1:3]

# can have a peek at the data -> boxplot(counts(deseq_dataset, normalized=TRUE))

counts(deseq_dataset)[1:6, 1:3]
#RNA counts are not normally distributed
#PCA is a parametric technique and assumes normality
#RNA counts are roughly log normally distributed (zero values may produce issues)
#can do a regularized log transformation OR variant stabilizing transformation
#they differ in speed (not so much functionally) - variant is faster

vst = varianceStabilizingTransformation(deseq_dataset)
# boxplot(assay(vst)) ---- need this 'assay' to extract the data from vst

plotPCA(vst, intgroup = 'conditions') + 
  theme_bw()


#---- importante
#PCA- the closer the two plots are together the more similar they are
#to one another
#PC explain a proportion of the variance in the data
#analyzes top 500 most variable genes (this algo)
#usually the first PC would explain around 50 % not 98 %
#we can see that the split here is separated by the different cell lines
#implies that they were too different to compare in the first place


#---- we need to split the deseq data table in two halves
#and the re-normalize

dds1 = deseq_dataset[,1:6]
dds2 = deseq_dataset[,7:12]

#last time we generated normalization the counts were not
#overwritten
#we use normalized=TRUE to look at it

dds1 = estimateSizeFactors(dds1)
dds2 = estimateSizeFactors(dds2)


counts(dds1, normalized=TRUE)[1:6,]

#re-do PCA
vst1 = varianceStabilizingTransformation(dds1)
plotPCA(vst1, intgroup='conditions')

#and for the other cell line 
vst2 = varianceStabilizingTransformation(dds2)
plotPCA(vst2, intgroup='conditions')

#now it looks good 
#dds1:first axis - variance because of gene expression
#mock samples are more varaible than infected samples 
#dds2: infected 

#sample distance
#PCA is not very quantitative (semi)
#you could use: clustering (supervised and un)
#R is good for clustering - statistical
#use base libraries

# can use heirarchal clustering
#vst object is in a datatype - summarized experiment
assay(vst1) # assay() to extract the data
dist(t(assay(vst1))) # needs a distance matrix use dist
# needs to be transposed first
h =hclust(dist(t(assay(vst1)))) 

plot(h) # will produce a clustered dendrogram

# can us k-means clustering for a supervised clustering

#k = kmeans(t(assay(vst1)), centers=2) # takes numeric data and number of clusters

#k$cluster # there you go

#there are higher dimension methods tisny, umap
#for dimensionality reduction 


# day 32 determining differentially expressed genes----

#all of the analyses up to this point showed good quality

# two more DESeq steps 

#2)assess variance
#3)pull out gene data

#DESeq
#1) estimate the size factors (normalization)
#2) estimate dispersion (sort of like variance)
#biol variance (highly expressed genes) + technical variance (low expressed genes)
#measure a sample of population instead of the whole pop
#in DESeq uses variance for 60k genes and submit it to shrinkage
#empirical Bayesian shrinkage
#3) apply statistics (Wald Test)



#dds1 = estimateDispersions(dds1)
#dds1

# sample_files=paste0(pull(sample_table, `Sample Name`), '/quant.sf') 
# 
# #give it names (otherwise it would be just numbered)
# 
# names(sample_files) = pull (sample_table, `Sample Name`)

#the 1:6 slice after import causes problems as some metadata still expects 4 conditions

#need to reimport sliced

sample_files_nhbe = sample_files[1:6]

count_data_nhbe = tximport(files = sample_files_nhbe,
                     type = 'salmon',
                     tx2gene = gene_map,
                     ignoreTxVersion = TRUE)

sample_table_nhbe = sample_table[1:6,]
sample_table_nhbe$conditions = factor(rep(c('mock', 'infected'), each=3), levels=c('mock', 'infected'))

# levels to have 'mock' first as DESEQ wants that 

dds_nhbe = DESeqDataSetFromTximport(txi = count_data_nhbe,
                                         colData= sample_table_nhbe,
                                         design =~conditions)

dds_nhbe = estimateSizeFactors(dds_nhbe)


# back to business now

dds_nhbe = estimateDispersions(dds_nhbe)

# have a look at them using DESEQ

plotDispEsts(dds_nhbe)

# as mean increases variance goes down in RNA-sEQ
# black points are gene wise estimations of dispersion
# blue points are the final dispersion estiamtes 'final'
# looks how we'd expect it to look

#Step 3 ----
dds_nhbe = nbinomWaldTest(dds_nhbe)

#Cheat 123 ---- dds_nhbe = DESeq(dds_nhbe) applies all three functions automatically

results_table=results(dds_nhbe) # pairwise comparison, no need to define anything

summary(results_table) # get a summary of the data

#DESEQ tries to preserve the most differently expressed genes
# low mean count filter will make you lose positive results

#LFC > 0 (up)       : 297, 1%   - in infected
#LFC < 0 (down)     : 168, 0.59% - in infected 

head(results_table)
View(as.data.frame(results_table))

# Day 33 ----
#results_table is a DataFrame and not data.frame
#let's make a copy as a data.frame - for functions

results_df = as.data.frame(results_table)

#check a gene to see if it has an outlier (for cooks distance)

plotCounts(dds_nhbe, gene = 'ENSG00000265794',
           intgroup='conditions')
#look at the outlier

#at the moment not all rows are informative of differential expression
#let's remove the filtered genes from our statistics

#any filtered gene has N/A
sum(complete.cases(results_df)) # sum() will sum up the TRUEs

#what follows is essentially a .dropna
filter_df1 = results_df[complete.cases(results_df),] # filter the columns 



#what filters do we want to apply?
#padj < 0.05 
# log2foldchange >1<-1


filter_df1$padj < 0.05

filter_df2 = filter_df1[filter_df1$padj < 0.05,]# comma because only filtering the rows and not columns

abs(filter_df2$log2FoldChange)  # gives absolute value (modulo)

abs(filter_df2$log2FoldChange) > 1

#all genes tested by DESeq
# have padj 0.05
#and have log2foldchange larger that 1 absolute
filter_df3 = filter_df2[abs(filter_df2$log2FoldChange) > 1,]
dim(filter_df3) # look at the shape, the rows shows us the number of
#differentially expressed genes

#day 34 ----
#TODO make 2-3 plots
#MA plot
plotMA(results_table)
#produce comaprison of mean of normalized counts and change
#genes with low count have high fold changes
#low count genes have high log changes but low p values
#technical variance
# colored genes are significantly differential genes

#Volcano plot
# X-axis log2FC, Y axis - p-value (log transformed -log10, the smaller p value the higher after transform)
#p-value tends to increase with log2FC
#we need 
#use filter_df1 (NAN dropped) contains all genes and adjuster p-value


#build it using GG plot

#make two filers on unfiltered data from yesterday to color the volcano plot 
filter_df1$test = filter_df1$padj < 0.05 & abs(filter_df1$log2FoldChange) > 1
# & combines conditions 
#above we added a 'test' column that says true or false based on matching the filter


#to label dots need to put row names into columns 
filter_df1 = rownames_to_column(filter_df1, var='ensgene')
#ENSG ids are not very human readable though 
#We will use biomart to achieve this
#can use the API install('biomaRt')

g = ggplot(filter_df1, aes(x=log2FoldChange, y=-log10(padj),
                           name=ensgene))+
  geom_point(aes(colour=test), size=1, alpha=0.3)+
  scale_colour_manual(values=c('black', 'red'))+
  geom_vline(xintercept=1, colour='green', linetype=3)+
  geom_vline(xintercept=-1, colour='green',linetype=3)+
  geom_hline(yintercept=-log10(0.05), colour='blue',linetype=3)+
  xlim(-3,3)+
  ylim(0,10)+
  theme_bw() + 
  theme(legend.position = 'none')
  
#make interactive plots by putting our 
ggplotly(g)


#bioMart----
listMarts()
#ensemble 100 uses gencode 34, we have 33 (ensemble 99)
ensembl99 = useEnsembl(biomart="ensembl", version = 99 ) #connect
View(listDatasets(ensembl99)) # 
ensembl99 = useDataset("hsapiens_gene_ensembl",
                       mart=ensembl99)
#we have just selected the human dataset
#what attributes do we have?
View(listAttributes(ensembl99))
View(listFilters(ensembl99))
#need to tick the boxes for attributes we want
getBM(attributes = c('ensembl_gene_id', 'ensembl_gene_id_version',
                     'ensembl_transcript_id',	
                     'ensembl_transcript_id_version',
                     'external_gene_name'),
                     filters = c('ensembl_gene_id'),
                     values = filter_df1$ensgene[1:6],
                     mart= ensembl99) # filter based on e.g. ENSG00000000000000003

#now let's get them all with useful to know attributes 
annotation = getBM(attributes = c('ensembl_gene_id', 'chromosome_name',
                                  'start_position',
                                  'end_position',
                                  'strand',
                                  'gene_biotype',
                                  'external_gene_name',
                                  'description'),
                   filters = c('ensembl_gene_id'),
                   values = filter_df1$ensgene,
                   mart= ensembl99)

View(annotation)

# let's merge this to our results
#use by=c() to say which columns are supposed to be the same
annotated_df = left_join(filter_df1, annotation, by=c('ensgene'='ensembl_gene_id') )

dim(annotated_df)
#15 in total

g = ggplot(annotated_df, aes(x=log2FoldChange, y=-log10(padj),
                           name=external_gene_name))+
  geom_point(aes(colour=test), size=1, alpha=0.3)+
  scale_colour_manual(values=c('black', 'red'))+
  geom_vline(xintercept=1, colour='green', linetype=3)+
  geom_vline(xintercept=-1, colour='green',linetype=3)+
  geom_hline(yintercept=-log10(0.05), colour='blue',linetype=3)+
  xlim(-3,3)+
  ylim(0,10)+
  theme_bw() + 
  theme(legend.position = 'none')

ggplotly(g)


#day 36 ----
#make a heatmap of differentially expressed genes


View(filter_df3)


dim(filter_df3)

annotated_df2 = annotated_df[annotated_df$padj < 0.05,]# comma because only filtering the rows and not columns

annotated_df3 = annotated_df2[abs(annotated_df2$log2FoldChange) > 1,]

View(annotated_df3)

degs = annotated_df3$ensgene

#get transformed counts 
vst_nhbe = varianceStabilizingTransformation(dds_nhbe)
#extract the vst counts
vst_nhbe_mat = assay(vst_nhbe)

#filter the data frame
data_for_hm = vst_nhbe_mat[degs,]
rownames(data_for_hm) = annotated_df3$external_gene_name

# get the heatmap
heatmap(data_for_hm) #scaling by row applied Z-score


#to pick a nice colour scheme you can use colorbrewer2.org

#display.brewer.all()

#look at this 
#extrapolate from 9 colours to 100 shades
greys = colorRampPalette(brewer.pal(9, name = "Greys"))(100)


pheatmap(data_for_hm, fontsize_row = 4, scale='row', color=greys,
         cutree_rows = 2,
         cutree_cols = 2) 

last_scheme = colorRampPalette(brewer.pal(7, "P"))



  #----Day37 functional analysis, gene ontology
#today we will do gene ontology analysis to understand what these tools do
#GO comprehensive dictionary of term related to biological systems
# pathways, functions, proteases, cell localization
# What we have now is a list of differentially expressed genes
#we now want to add a list of functions to this list

#we will use hypergeometric test to see if there is statistical enrichment 
#of a specific proces

#visualize using 'pot of balls'
# need to install clusterProfiler from BiocManager
# and BiocManager::install('org.Hs.eg.db')


# use enrichGO command
# we are going to gene gene names from biomart
#
ego = enrichGO(
  gene = ,
  OrgDb = ,
  ont = "MF",
  universe = 
)

# just nete the entrez gene ID
ent_gene = getBM(attributes = c('entrezgene_id'),
                              filters = c('ensembl_gene_id'),
                              values = annotated_df3$ensgene,
                              mart= ensembl99)
ent_gene
#lets make it into a vector
ent_gene= ent_gene$entrezgene_id
ent_gene= as.character(ent_gene)


#we also need to define the pool to which we compare the dataset to
#everything that passed before testing for differential gene expression
ent_uni = getBM(attributes = c('entrezgene_id'),
                 filters = c('ensembl_gene_id'),
                 values = annotated_df$ensgene,
                 mart= ensembl99)
ent_uni = ent_uni$entrezgene_id
ent_uni = as.character(ent_uni)

ego = enrichGO(
  gene = ent_gene,
  OrgDb = org.Hs.eg.db ,
  ont = "BP",
  universe = ent_uni
)
#it's referencing all of the gene ontology data 
#enrichGO results can be summarized as a table

summary(ego)

View(as.data.frame(ego))

#quick visualization of top GO terms
barplot(ego)
barplot(ego, showCategory=20)
dotplot(ego)


ekg = enrichKEGG (gene = ent_gene,
                  universe = ent_uni)

#look at pathway enrichment 
View(as.data.frame(ekg))

write_tsv(annotated_df3, "filtered_nhbe_results.txt")



