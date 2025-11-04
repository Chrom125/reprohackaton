
#!/usr/bin/env Rscript

############################## Arguments Parsing ################################
library(argparse)

#Defining the arguments parser
parser = ArgumentParser(description= 'This script performs differential expression analysis using DESeq2 package 
and generates MA plots for the entire dataset and translation-related genes')

# Defining the arguments
parser$add_argument('--countTable', '-c', help= 'Path to the count table file', required= TRUE)
parser$add_argument('--metadata', '-m', help= 'Path to the metadata file (tsv format)', required= TRUE)
parser$add_argument('--geneNames', '-g', help= 'Path to the gene names file (Excel format)', required= TRUE)
parser$add_argument('--outputDir', '-o', help= 'Path to the output directory', required= FALSE, default= './')

# Parsing the arguments
xargs = parser$parse_args()

############################## D.E Analysis ####################################
library(DESeq2)
###### Step 1 : Data preparation

#-------------------------------------------------------------------------------
# Count table
#-------------------------------------------------------------------------------
count_table = read.table(file = xargs$countTable, header = T, sep = "\t", 
                         row.names = 1)
head(count_table) 

#Total number of genes
number_of_genes = dim(count_table)[1]
number_of_genes
#2967

#-------------------------------------------------------------------------------
# Metadata table
#-------------------------------------------------------------------------------
metadata = read.table(xargs$metadata, header = T, sep ="\t")
rownames(metadata) = metadata$Sample
#It's important to have rownames(metadata) in colnames(count_data)
#It helps DESeq2 to match samples with their description in the metadata table

#-------------------------------------------------------------------------------
# Count table in the DESeq2 format
#-------------------------------------------------------------------------------

dds = DESeqDataSetFromMatrix(countData = as.matrix(count_table),
                             colData = metadata,
                             design = ~Label)

#Defining the reference for evaluating the differential expression
dds$Label = relevel(dds$Label, "Control")



##### Step 2 : Running D.E Analysis
dds.DE = DESeq(dds)

#Accessing the normalized counts
normalized_counts <- counts(dds.DE, normalized = TRUE)
normalized_counts



#### Step3 : Results analysis
DE.results = results(dds.DE)
DE.results

# log2 fold change (MLE): Label Persister vs Control 
# Wald test p-value: Label Persister vs Control 
# DataFrame with 2967 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   SAOUHSC_00001 13042.560       0.274100  0.275688  0.994238 3.20107e-01 4.13735e-01
# SAOUHSC_00002  9610.983       0.448343  0.180917  2.478164 1.32060e-02 2.67522e-02
# SAOUHSC_00003   822.104       0.663142  0.359727  1.843457 6.52623e-02 1.07973e-01
# SAOUHSC_00004  6998.784       1.309560  0.280318  4.671698 2.98719e-06 1.40098e-05
# SAOUHSC_00005 16195.254       1.288328  0.333230  3.866187 1.10550e-04 3.68011e-04


#### DE analysis details for differentially expressed genes only
padj.treshold = 0.05
differentially_expressed_genes = DE.results[which(DE.results$padj < padj.treshold),]

#### Number of differentially expressed of genes

dim(differentially_expressed_genes)[1]
#1483

#### Number of genes up-regulated
dim(differentially_expressed_genes[differentially_expressed_genes$log2FoldChange > 0,])[1]
#712

#### Number of genes down-regulated
dim(differentially_expressed_genes[differentially_expressed_genes$log2FoldChange < 0,])[1]
#771



##### MA plot for the entire RNAseq DataSet
library(ggplot2)
library(scales) # for formatting axis graduations

#Collecting BaseMean, lFC and padj for all genes
dataf.DE.results = as.data.frame(DE.results[,c("baseMean","log2FoldChange","padj")])
#Adding GeneID columns to dataf.DE.results
dataf.DE.results = data.frame(GeneID = rownames(dataf.DE.results), dataf.DE.results)

#Adding info on the significance of differential expression
dataf.DE.results$signif = ifelse(!is.na(dataf.DE.results$padj) & dataf.DE.results$padj < padj.treshold,
  "Significant",               # significantly differentially expressed
  "Non-Significant"            # Not significantly differentially expressed
)

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



##### MA plot for translation genes


#-------------------------------------------------------------------------------
# Geneid geneName equivalence table (from aureowiki)
#-------------------------------------------------------------------------------
library(readxl)
geneName_geneID_eq = read_excel(xargs$geneNames)
#Renaming colums
colnames(geneName_geneID_eq) = c("GeneID", "Organism", "GeneName", "Product")

#number of lines
n = dim(geneName_geneID_eq)[1]

#Removing empty the last three  rows which do not contain gene info
geneName_geneID_eq = geneName_geneID_eq[1:(n-3),]

#Keeping only genes feature type:  with GeneID like "SAOUHSC_*" 
geneName_geneID_eq = geneName_geneID_eq[grepl("^SAOUHSC_", geneName_geneID_eq$GeneID), ]

# Rapid check of gene content and differences with our count table
setdiff(geneName_geneID_eq$GeneID, rownames(count_table))

# annotations present in  in Aureowiki and absent in our table
# [1] "SAOUHSC_00381a"              "SAOUHSC_00411.1"             "SAOUHSC_00411.2"            
# [4] "SAOUHSC_00411.3"             "SAOUHSC_00411.4"             "SAOUHSC_00630.1"            
# [7] "SAOUHSC_01055a"              "SAOUHSC_01139/SAOUHSC_01141" "SAOUHSC_1307a"              
# [10] "SAOUHSC_1342a"               "SAOUHSC_01761a"              "SAOUHSC_T00046"             
# [13] "SAOUHSC_T00047"              "SAOUHSC_02512a"              "SAOUHSC_03037a"             
# [16] "SAOUHSC_3042a"              

# annotations present in our table and absent in Aureowiki
setdiff(rownames(count_table),geneName_geneID_eq$GeneID)
# [1] "SAOUHSC_00380"  "SAOUHSC_00595"  "SAOUHSC_00630"  "SAOUHSC_00974"  "SAOUHSC_01057"  "SAOUHSC_01139" 
# [7] "SAOUHSC_01140"  "SAOUHSC_01141"  "SAOUHSC_01302"  "SAOUHSC_01309"  "SAOUHSC_01342"  "SAOUHSC_01344" 
# [13] "SAOUHSC_01410"  "SAOUHSC_01648"  "SAOUHSC_01725"  "SAOUHSC_01804"  "SAOUHSC_01881"  "SAOUHSC_01906" 
# [19] "SAOUHSC_A01910" "SAOUHSC_02231"  "SAOUHSC_02302"  "SAOUHSC_02438"  "SAOUHSC_02440"  "SAOUHSC_02473" 
# [25] "SAOUHSC_02556"  "SAOUHSC_02633"  "SAOUHSC_02637"  "SAOUHSC_02662"  "SAOUHSC_02673"  "SAOUHSC_02787" 
# [31] "SAOUHSC_02814"  "SAOUHSC_02903"  "SAOUHSC_03038"  "SAOUHSC_03043" 


#Specific verification for the "SAOUHSC_01139/SAOUHSC_01141" entry

#In our table
dataf.DE.results[grepl("^SAOUHSC_01139", rownames(dataf.DE.results)), ]
dataf.DE.results[grepl("^SAOUHSC_01141", rownames(dataf.DE.results)), ]

#In the Aureowiki table
geneName_geneID_eq[grepl("^SAOUHSC_01139/SAOUHSC_01141", geneName_geneID_eq$GeneID), ]
# # A tibble: 1 × 4
# GeneID                      Organism           GeneName Product                      
# <chr>                       <chr>              <chr>    <chr>                        
#   1 SAOUHSC_01139/SAOUHSC_01141 S. aureus NCTC8325 bshC     Putative cysteine ligase BshC


#Adding geneName info to the DESeq2 results table (dataf.DE.results)
dataf.DE.results.extended = merge(dataf.DE.results, geneName_geneID_eq, by = "GeneID", all.x = TRUE)
dataf.DE.results.extended$GeneName[dataf.DE.results.extended$GeneID=="SAOUHSC_01139"] = "bshC"
dataf.DE.results.extended$GeneName[dataf.DE.results.extended$GeneID=="SAOUHSC_01141"] = "bshC"


#-------------------------------------------------------------------------------
# Metabolic pathways for S.aureus genes (from KEGG)
#-------------------------------------------------------------------------------
library(KEGGREST)

#List of KEGG Database
listDatabases()

#Correspondance table between S. aureus (code sao) genes and their KEGG BRITE functional hierarchy ID
GeneID_BRITE.ID= keggLink(target = "brite", source = "sao") # vector of KEGG BRITE functional hierarchy (with the genes as names)

GeneID = sub(pattern = "^sao:",  replacement = "", names(GeneID_BRITE.ID)) #genes ID without the prefix "sao:"

BRITE.ID = sub(pattern = "^br:",  replacement = "", GeneID_BRITE.ID) #KEGG BRITE functional hierarchy ID without the prefix "br:"

GeneID_BRITE.ID= data.frame(GeneID, BRITE.ID) # table GeneID | KEGG BRITE ID


#Aggregating all the hierarchies for each gene, since a gene may be active in many functional pathway
GeneID_BRITE.ID = aggregate(. ~ GeneID, data = GeneID_BRITE.ID,
                             FUN = function(x) paste(unique(x), collapse = ";"))

#Adding KEGG BRITE functional hierarchy ID info to the DESeq2 results table (dataf.DE.results.extended)

dataf.DE.results.extended = merge(dataf.DE.results.extended, GeneID_BRITE.ID,
                                  by= "GeneID", all.x = TRUE)

## MA plot drawing for translation related genes
#===========> In KEGG database, S. aureus translation genes are grouped in 03 BRITE functional hierarchies
#===========>sao03011  Ribosome
#===========>sao03009  Ribosome biogenesis
#===========>sao03016  Transfer RNA biogenesis
#===========>sao03012  Translation factors


#filter for all translation related genes
filter_translation_genes = !is.na(dataf.DE.results.extended$BRITE.ID) &
  (grepl("(^|;)sao03011(;|$)", dataf.DE.results.extended$BRITE.ID) |
   grepl("(^|;)sao03009(;|$)", dataf.DE.results.extended$BRITE.ID) |
   grepl("(^|;)sao03016(;|$)", dataf.DE.results.extended$BRITE.ID) |
   grepl("(^|;)sao03012(;|$)", dataf.DE.results.extended$BRITE.ID))

#filter for AA_tRNA_synthetases only
filter_AA_tRNA_synthetases = !is.na(dataf.DE.results.extended$BRITE.ID) &
  grepl("(^|;)sao03016(;|$)", dataf.DE.results.extended$BRITE.ID) &
  grepl("-tRNA synthetase", dataf.DE.results.extended$Product)

#Typical genes in translation  
typical_members = subset(
  dataf.DE.results.extended,
  GeneName %in% c("pth", "infA", "infB", "infC", "frr", "tsf")
)

# Offsets for plot annotation
typical_members$dx <- c(-1.8, 0.8, 0.6, 2.1, 2, 2.5)
typical_members$dy <- c(-0.8, 1.3, 1.7, 0, 0, -0.7)

f2 = ggplot(data = dataf.DE.results.extended[filter_translation_genes,], 
            aes(x = log(baseMean, base = 2) , y = log2FoldChange, color = signif))+
     geom_point(alpha = 0.8,size = 0.9)+ #Plotting all translation related genes
     geom_point(data = dataf.DE.results.extended[filter_AA_tRNA_synthetases, ],
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
