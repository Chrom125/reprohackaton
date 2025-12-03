#!/usr/bin/env Rscript

library(DESeq2) # For differential expression analysis

##################### Step 1 : Data preparation ################################
################################################################################

##### Count table
count_table = read.delim(file = "results/comparison_with_authors/authors_counts.xls", stringsAsFactors=TRUE)
count_table = count_table[1:2967,c(2,6:11)] #Selecting only gene counts columns
colnames(count_table)[1] = "GeneID" #Renaming first column to GeneID
rownames(count_table) = count_table$GeneID #Setting GeneID as rownames
count_table$GeneID = NULL #Removing GeneID column as it's now in rownames
#Total number of genes
number_of_genes = dim(count_table)[1]
cat("Number of genes in the count table : ", number_of_genes, "\n")

##### Metadata table

#Specifying the sample nature via a variable called Label
metadata = data.frame(Label = factor(c(rep("Control",3),rep("Persister", 3))))
#It's important to have rownames(metadata) in colnames(count_data)
#It helps DESeq2 to match samples with their description in the metadata table
rownames(metadata) = colnames(count_table)

##### Count table in DESeq2 Format

dds = DESeqDataSetFromMatrix(countData = as.matrix(count_table),
                             colData = metadata,
                             design = ~Label)

#Defining the reference for evaluating the differential expression
dds$Label = relevel(dds$Label, "Control")



################################################################################
##################### Step 3 : Running DESeq analysis ##########################
################################################################################
dds.DE = DESeq(dds)


################################################################################
######################## Step 4 : Results analysis #############################
################################################################################

################################## Preparing results table #####################

# Adjusted p-value cutoff to consider a gene as differentially expressed
padj.cutoff = 0.05

##### DEseq2 results extraction
DE.results = results(dds.DE, alpha = padj.cutoff)
cat("log2 fold change (MLE): Label Persister vs Control\n")
cat("Wald test p-value: Label Persister vs Control\n")
#Results summary (Number of up/down regulated genes)
summary(DE.results)


#### Collecting BaseMean, lFC and padj for all genes
dataf.DE.results = as.data.frame(DE.results[,c("baseMean","log2FoldChange","padj")])
#Adding GeneID columns to dataf.DE.results
dataf.DE.results = data.frame(GeneID = rownames(dataf.DE.results), dataf.DE.results)


#### Writing DESeq2 results table for all genes
write.table(dataf.DE.results,
            file = file.path("results/comparison_with_authors/DESeq2_results_authors.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
