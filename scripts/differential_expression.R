
#!/usr/bin/env Rscript
############################## D.E Analysis ####################################
################################################################################
#Useful libraries
library(ggplot2) # For data visualization
library(ggrepel)  # For better text label placement in ggplot2
library(DESeq2) # For differential expression analysis
library(readxl) # For reading Excel files
library(scales) # For formatting axis graduations

##################### Step 1 : Data preparation ################################
################################################################################

##### Count table
count_table = read.table(file = "results/featurecounts/processed_counts.txt", header = T, sep = "\t", 
                         row.names = 1)
                           
#Total number of genes
number_of_genes = dim(count_table)[1]
cat("Number of genes in the count table : ", number_of_genes, "\n")

##### Metadata table

#Specifying the sample nature via a variable called Label
metadata = data.frame(Label = factor(c(rep("Persister",3),rep("Control", 3))))
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
##################### Step 2 : Clustering using PCA ############################
################################################################################

#### PCA plot

#Applying a variance stabilizing transformation to raw counts data
transformed_counts = vst(dds, blind = FALSE)

f1 = plotPCA(transformed_counts, intgroup="Label") +
     labs(color = "Sample",
          title = "Principal component analysis")
f1 = f1 + geom_text_repel(aes(label = name), size = 3, nudge_x = 0.5,         
                          nudge_y = 0.5,        
                          force = 2,
                          segment.size = 0.2 )
ggsave("results/DESeq2_Analysis/PCA_plot.png", plot = f1, width = 6, height = 5, dpi = 600)


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

#### Adding info on the significance of differential expression
dataf.DE.results$signif = ifelse(!is.na(dataf.DE.results$padj) & dataf.DE.results$padj < padj.cutoff,
                                  "Significant",     # significantly differentially expressed
                                 "Non-Significant"   # Not significantly differentially expressed
)

#### Adding gene name information from Aureowiki
cat("About gene name information from Aureowiki: S. aureus Strain NCTC8325\n")
geneName_geneID_eq = read_excel("Data/GeneSpecificInformation_NCTC8325.xlsx")
#Renaming colums
colnames(geneName_geneID_eq) = c("GeneID", "Organism", "GeneName", "Product")
#number of lines
n = dim(geneName_geneID_eq)[1]
#Total number of features annotated in Aureowiki
cat("Total number of features annotated in aureowiki", n, "\n")

#Removing  the last three empty rows which do not contain gene info
geneName_geneID_eq = geneName_geneID_eq[1:(n-3),]

#Keeping only genes feature type:  with GeneID like "SAOUHSC_*" 
geneName_geneID_eq = geneName_geneID_eq[grepl("^SAOUHSC_", geneName_geneID_eq$GeneID), ]

#Total number of genes features annotated in Aureowiki
cat("Total number of gene features annotated in aureowiki (with GeneID like SAOUHSC_*)  :", dim(geneName_geneID_eq)[1], "\n\n")

##Rapid check of gene content and differences with our count table
# Gene annotations present in Aureowiki and absent in our count table
cat("Annotations present in Aureowiki and absent in our count table :\n\n")
setdiff(geneName_geneID_eq$GeneID, rownames(count_table))
cat("\n\n")
# Gene annotations present in our table and absent in Aureowiki
cat("Annotations present in our count table and absent in Aureowiki :\n\n")
setdiff(rownames(count_table),geneName_geneID_eq$GeneID)
cat("\n\n")
# Specific verification for the "SAOUHSC_01139/SAOUHSC_01141" entry
cat("The entry 'SAOUHSC_01139/SAOUHSC_01141' in Aureowiki corresponds to two separate genes in our count table:\n")

#In the Aureowiki table
cat("In Aureowiki table:\n")
as.data.frame(geneName_geneID_eq[grepl("^SAOUHSC_01139/SAOUHSC_01141", geneName_geneID_eq$GeneID), ])
cat("\n\n")

#In our table
cat("In our count table:\n")
count_table[grepl("^SAOUHSC_01139", rownames(count_table)), ]
cat("\n")
count_table[grepl("^SAOUHSC_01141", rownames(count_table)), ]
cat("\n")

#### Adding geneName info to the DESeq2 results table (dataf.DE.results)
dataf.DE.results = merge(dataf.DE.results, geneName_geneID_eq, by = "GeneID", all.x = TRUE)
dataf.DE.results$GeneName[dataf.DE.results$GeneID=="SAOUHSC_01139"] = "bshC"
dataf.DE.results$GeneName[dataf.DE.results$GeneID=="SAOUHSC_01141"] = "bshC"


#KEGG BRITE Gene functional annotation file
GeneID_BRITE.ID = read.delim("results/Gene_KEGG_functional_annotation/KEGG_BRITE_functional_annotation.tsv")

#Adding KEGG BRITE functional hierarchy ID to the DESeq2 results table (dataf.DE.results)
dataf.DE.results = merge(dataf.DE.results, GeneID_BRITE.ID, by= "GeneID", all.x = TRUE)
cat("\n About S. aureus strain NCTC8325 gene functional annotation in KEGG BRITE: \n")
cat("Total number of genes with KEGG BRITE functional hierarchy annotation : ", dim(GeneID_BRITE.ID)[1], "\n")

#### Writing DESeq2 results table for all genes
write.table(dataf.DE.results, file = "results/DESeq2_Analysis/DESeq2_results_all_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

