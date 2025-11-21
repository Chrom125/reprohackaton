
#!/usr/bin/env Rscript

library(argparse) # For parsing command line arguments
############################## Arguments Parsing ################################


#Defining the arguments parser
parser = ArgumentParser(description= 'This script performs differential expression 
analysis using DESeq2 package 
and generates MA plots for the entire dataset and translation-related genes')

# Defining the arguments
parser$add_argument('--countTable', '-c', help= 'Path to the count table file', 
                    required= TRUE)
parser$add_argument('--geneNames', '-g', help= 'Path to the gene names file (Excel format)', 
                    required= TRUE)

parser$add_argument('--geneFunction', '-gF', help= 'Path to the gene KEGG functionnal annotation file', 
                    required= TRUE)

parser$add_argument('--outputDir', '-o', help= 'Path to the output directory', 
                    required= FALSE, default= './')

# Parsing the arguments
if (interactive()) {
    # Valeurs par défaut pour travail interactif (adaptez les chemins)
    xargs = list(
      countTable = "results/featurecounts/processed_counts.txt",
      geneNames  = "GeneSpecificInformation_NCTC8325.xlsx",
      geneFunction = "KEGG_BRITE_functional_annotation.tsv",
      outputDir  = "./results/DESeq2_bowtie2/"
    )
} else {
    xargs = parser$parse_args()
}

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
count_table = read.table(file = xargs$countTable, header = T, sep = "\t", 
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
ggsave(file.path(xargs$outputDir,"PCA_plot.png"), plot = f1, width = 6, height = 5, dpi = 600)


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
geneName_geneID_eq = read_excel(xargs$geneNames)
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
GeneID_BRITE.ID = read.delim(xargs$geneFunction)

#Adding KEGG BRITE functional hierarchy ID to the DESeq2 results table (dataf.DE.results)
dataf.DE.results = merge(dataf.DE.results, GeneID_BRITE.ID,
                                  by= "GeneID", all.x = TRUE)
cat("\n About S. aureus strain NCTC8325 gene functional annotation in KEGG BRITE: \n")
cat("Total number of genes with KEGG BRITE functional hierarchy annotation : ",
    dim(GeneID_BRITE.ID)[1], "\n")

#### Writing DESeq2 results table for all genes
write.table(dataf.DE.results,
            file = file.path(xargs$outputDir, "DESeq2_results_all_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)



##################################### Plots ###############################

##### MA plot for the entire RNAseq DataSet

#Handling out of range values
limit = 4 # max absolute value for y (log2FoldChange) axis
offset = 0.2 # the amount to offset the points from the limit.
dataf.DE.results$lfc2 = ifelse(abs(dataf.DE.results$log2FoldChange) > limit, sign(dataf.DE.results$log2FoldChange)*limit + sign(dataf.DE.results$log2FoldChange)*offset, 
                               dataf.DE.results$log2FoldChange) # cap x-values at -lim or lim depending on which limit they exceed
dataf.DE.results$flag = abs(dataf.DE.results$log2FoldChange) > limit # add a flag variable to help filter for points outside the limits

f1 = ggplot(dataf.DE.results, aes(x = baseMean , y = lfc2, color = signif))+
     geom_point(data = subset(dataf.DE.results, !flag), 
                alpha = 0.6,size = 0.8)+ # points inside the limits
     geom_point(data = subset(dataf.DE.results, flag==TRUE & log2FoldChange < -1*limit), 
                alpha = 0.6,size = 0.8, shape = 25, 
                aes(fill = signif))+ # points outside the limits with negative lfc
     geom_point(data = subset(dataf.DE.results, flag==TRUE & log2FoldChange > limit), 
                alpha = 0.6,size = 0.8, shape = 24, 
                aes(fill = signif))+ # points outside the limits with positive lfc
     scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + #logarithmic scale on x-axis                                    
     scale_color_manual(values = c("Non-Significant" = "black",
                                "Significant" = "red")) +
     geom_hline(yintercept = 0, linetype = "dashed") +
     coord_cartesian(ylim = c(-4, 4))+
     labs(
          x = "Mean of normalized counts",
          y = expression(log[2]~"fold change"))+theme(legend.position = "none")

# Saving
ggsave(file.path(xargs$outputDir,"MA_plot_all_genes.png"), plot = f1, width = 6, height = 5, dpi = 300)


######## Volcano plot

#Info on differential expression for all genes
dataf.DE.results$diffexpression = 'Not significant'
dataf.DE.results$diffexpression[dataf.DE.results$log2FoldChange>0 & dataf.DE.results$signif == "Significant"] = "Up-regulated"
dataf.DE.results$diffexpression[dataf.DE.results$log2FoldChange<0 & dataf.DE.results$signif == "Significant"] = "Down-regulated"

#Top 10 Genes with the most significant differential expression
top10degs = head(dataf.DE.results[order(dataf.DE.results$padj),],10)

volcano = ggplot(data = dataf.DE.results, aes(x = log2FoldChange, y = -log10(padj),
                                              col = diffexpression))+
  geom_vline(xintercept = c(-0.6, 0.6), col = 'gray', linetype = 'dashed')+
  geom_hline(yintercept = -log10(c(0.05)), col = 'gray', linetype = 'dashed')+
  geom_point()+
  scale_color_manual(values = c("blue","grey","red"))+ #color of the diffexpression value
  scale_x_continuous(breaks = seq(-7,7,1))+
  geom_label_repel(data = top10degs, aes(label = GeneName), show.legend = FALSE,
                   min.segment.length = 0)+
  labs(color = "Differential expression",
       x = expression(log[2]~"Fold Change"),
       y = expression("-log"[10]~"(adjusted p-value)"))+ 
  theme_classic()
ggsave(file.path(xargs$outputDir,"Volcanoplot.png"), plot = volcano, dpi = 600)

##### MA plot for translation genes

#===========> In KEGG database, S. aureus translation genes are grouped in 03 BRITE functional hierarchies
#===========>sao03011  Ribosome
#===========>sao03009  Ribosome biogenesis
#===========>sao03016  Transfer RNA biogenesis
#===========>sao03012  Translation factors


#filter for all translation related genes
filter_translation_genes = !is.na(dataf.DE.results$BRITE.ID) &
  (grepl("(^|;)sao03011(;|$)", dataf.DE.results$BRITE.ID) |
   grepl("(^|;)sao03009(;|$)", dataf.DE.results$BRITE.ID) |
   grepl("(^|;)sao03016(;|$)", dataf.DE.results$BRITE.ID) |
   grepl("(^|;)sao03012(;|$)", dataf.DE.results$BRITE.ID))

#filter for AA_tRNA_synthetases only
filter_AA_tRNA_synthetases = !is.na(dataf.DE.results$BRITE.ID) &
  grepl("(^|;)sao03016(;|$)", dataf.DE.results$BRITE.ID) &
  grepl("-tRNA synthetase", dataf.DE.results$Product)

#Typical genes in translation  
typical_members = subset(
  dataf.DE.results,
  GeneName %in% c("pth", "infA", "infB", "infC", "frr", "tsf")
)

# Offsets for plot annotation
typical_members$dx <- c(-1.8, 0.8, 0.6, 2.1, 2, 2.5)
typical_members$dy <- c(-0.8, 1.3, 1.7, 0, 0, -0.7)

f2 = ggplot(data = dataf.DE.results[filter_translation_genes,], 
            aes(x = log(baseMean, base = 2) , y = log2FoldChange, color = signif))+
     geom_point(alpha = 0.8,size = 0.9)+ #Plotting all translation related genes
     geom_point(data = dataf.DE.results[filter_AA_tRNA_synthetases, ],
             aes(shape = "AA_tRNA_synthetases"), 
             color = "black", size = 1, stroke = 1.2) + #Plotting AA_tRNA_synthetases genes
    # color scale for significant and non-significant genes
     scale_color_manual(values = c("Non-Significant" = "grey",
                                "Significant" = "red")) +
    # shape scale for AA_tRNA_synthetases genes
     scale_shape_manual(
            values = c("AA_tRNA_synthetases" = 1))+
    # axis graduations
     scale_x_continuous(breaks = seq(0, 20, by = 2))+
     scale_y_continuous(breaks = seq(-6, 5, by = 1))+
    #Annotation segments for typical translation genes
     geom_segment(
      data = typical_members,
      aes(x = log(baseMean, base = 2)+0.05, y = log2FoldChange,
          xend = log(baseMean, base = 2) + dx, yend = log2FoldChange + dy),
          linewidth = 0.8, color = "black")+
  
     #Annotation text for typical translation genes
     geom_text(
      data = typical_members,
      aes(x = log(baseMean, base = 2)+ 1.35*dx, y = log2FoldChange + 1.15*dy, label = GeneName),
      inherit.aes = FALSE,
      fontface = "italic", size = 4) +
    # Horizontal line at y=0
     geom_hline(yintercept = 0, linetype = "dashed") +
  
     #Axis names + suppression of legend titles
     labs(
       x = expression(log[2]~"base Mean"),
       y = expression(log[2]~"fold change"), shape = NULL, color = NULL)+
  
     #Plot Graphic theme
     theme_classic() +
     theme(
      panel.border = element_rect(color = "black", linewidth = 2.5))

    
# Saving
ggsave(file.path(xargs$outputDir,"MA_plot_translation_genes.png"), plot = f2, width = 6, height = 5, dpi = 300)