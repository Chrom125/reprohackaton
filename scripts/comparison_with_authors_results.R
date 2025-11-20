################################################### Argument parsing #####################################################
library(argparse) # For argument parsing
#Defining the arguments parser
parser = ArgumentParser(description= "This script performs a comparison between authors's differential expression analysis results and ours by generating boxplots and an UpSet plot.")

# Defining the arguments
parser$add_argument('--authorsResults', '-aR', help= 'Path to authors DESeq2 results file', required= TRUE)
parser$add_argument('--reproducedResults', '-rR', help= 'Path to the reproduced DESeq2 results file', required= TRUE)
parser$add_argument('--outputDir', '-o', help= 'Path to the output directory', required= FALSE, default= './')

# Parsing the arguments
xargs = parser$parse_args()

############################################ Comparison with authors results #####################################################
#Libraries
library(scales) # Plot axis scales handling
library(ggplot2) # For Data Visualisation
library(ComplexHeatmap) #Library for Upset plot drawing
library(ragg) #Library to export plots
library(patchwork) # For plot organization

###### Collecting DESeq results for differentially expressed genes


###### Results from authors 

#Counts and DE analysis results from the authors
GSE139659_IPvsctrl.complete = read.delim(xargs$authorsResults, stringsAsFactors=TRUE)

dataf.DE.results.authors = GSE139659_IPvsctrl.complete[1:2967,c(2,12:17)]
colnames(dataf.DE.results.authors)[1] = "GeneID" 

#Differentially expressed genes
padj.cutoff = 0.05
differentially_expressed_genes.authors = dataf.DE.results.authors[which(dataf.DE.results.authors$padj < padj.cutoff),]

###### Results from us
#Counts and DE analysis results from us
dataf.DE.results.us = read.delim(xargs$reproducedResults, 
                                 stringsAsFactors=TRUE)

#Differentially expressed genes
differentially_expressed_genes.us = dataf.DE.results.us[which(dataf.DE.results.us$padj < padj.cutoff),]


######### Comparison of logFC, BaseMean for All genes between
######### us and authors

authors = dataf.DE.results.authors
#Column specifying that results are from authors
authors$From = rep("Authors", length(authors$GeneID))

us = dataf.DE.results.us
#Column specifying that results are from us
us$From = rep("Us", length(us$GeneID))

#Binding the two previous table
dataf.DE.results.all = rbind(authors[,c("GeneID", "baseMean","log2FoldChange","padj","From")],
                             us[,c("GeneID", "baseMean","log2FoldChange","padj", "From")])

#Adding Info on differentially expressed status
dataf.DE.results.all$Status = ifelse(dataf.DE.results.all$padj<padj.cutoff, 
                                     "Differentially Expressed", "Not Differentially Expressed")

#Boxplot of baseMean for all genes (Us vs Authors)

p1 = ggplot(dataf.DE.results.all,
                 aes(x = From, y = baseMean)) +
  geom_boxplot(outlier.size = 0.5, color = "black", fill = "pink") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  coord_flip() +
  labs(x = NULL,
       y = "Mean of normalized counts",
       title = "All genes (Us vs Authors)") +
  theme_bw()

#Boxplot of baseMean for differentially expressed (D.E)
#genes (Us vs Authors)

p2 = ggplot(subset(dataf.DE.results.all, Status == "Differentially Expressed"),
            aes(x = From, y = baseMean)) +
  geom_boxplot(outlier.size = 0.5, color = "black", fill = "grey") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  coord_flip() +
  labs(x = NULL,
       y = "Mean of normalized counts",
       title = "D.E genes (Us vs Authors) ") +
  theme_bw()


#Boxplot of logFC for all genes (Us vs Authors)

p3 = ggplot(dataf.DE.results.all,
            aes(x = From, y = log2FoldChange)) +
  geom_boxplot(outlier.size = 0.5, color = "black", fill = "pink")+
  coord_flip() +
  labs(x = NULL,
       y = expression(log[2]~"Fold Change"),
       title = "All genes (Us vs Authors)") +
  theme_bw()



#Boxplot of logFC for D.E (Us vs Authors)
p4 = ggplot(subset(dataf.DE.results.all, Status=="Differentially Expressed"),
            aes(x = From, y = log2FoldChange)) +
  geom_boxplot(outlier.size = 0.5, color = "black", fill = "grey") +
  coord_flip() +
  labs(x = NULL,
       y = expression(log[2]~"Fold Change"),
       title ="D.E genes (Us vs Authors)") +
  theme_bw()

pf = ((p1+p3)/(p2+p4)) + plot_annotation (tag_levels = "a", tag_suffix = ".")

ggsave(file.path(xargs$outputDir, "Boxplot_LFC_BaseMean_All_genes_and_DE_genes_Us_vs_Authors.png"),
       plot = pf, dpi = 600)


################### ANALYSIS of overlaps between our results and authors's

# Adding columns to the previous data frames indicating 
#if genes are up-regulated or down-regulated : 1 if True, 0 if not

#Up-regulated genes
differentially_expressed_genes.authors$up.regulated.authors = ifelse(differentially_expressed_genes.authors$log2FoldChange>0,
                                               1,0)
differentially_expressed_genes.us$up.regulated.us = ifelse(differentially_expressed_genes.us$log2FoldChange>0,1,0)

#Down-regulated genes
differentially_expressed_genes.authors$down.regulated.authors = ifelse(differentially_expressed_genes.authors$log2FoldChange<0,
                                               1,0)
differentially_expressed_genes.us$down.regulated.us = ifelse(differentially_expressed_genes.us$log2FoldChange<0,1,0)


#### Merging up and down regulated indications from authors and us
DE.genes.set = merge(differentially_expressed_genes.authors[,c("GeneID","up.regulated.authors","down.regulated.authors")], 
                     differentially_expressed_genes.us[,c("GeneID","up.regulated.us","down.regulated.us")], 
                     all = T, by="GeneID")

#Replacing all NA values by zeros in the previous gene set
DE.genes.set[is.na(DE.genes.set)] = 0


######### UPSET PLOT

#Combination matrix from genes_set
comb_matrix = make_comb_mat(DE.genes.set, mode = "distinct")

#Renaming genes set
set_name(comb_matrix) = c("Up-regulated (Authors)", "Down-regulated (Authors)",
                          "Up-regulated (Us)","Down-regulated (Us)")


f1= UpSet(comb_matrix,
          #Ordering the bars by decreasing intersection size
          comb_order = order(-comb_size(comb_matrix)),
          #Annotation of the barplot ath the right of the dot matrix
          right_annotation = rowAnnotation(
            "Gene sets size" = anno_barplot(set_size(comb_matrix), 
                                            ylim = c(0, max(set_size(comb_matrix))*1.1),
                                            add_numbers = TRUE,
                                            border = FALSE,
                                            gp = gpar(fill = "#5DA5DA", col = NA), 
                                            width = unit(6, "cm")), 
            annotation_name_side = "bottom", 
            annotation_name_rot = 0,
            annotation_name_gp = gpar(fontsize = 8),
            annotation_name_offset = unit(1.2, "cm")),
          
          #Annotation of the barplot ath the top of the dot matrix
          top_annotation = HeatmapAnnotation(
            "Intersection size\n(Number of genes)" = anno_barplot(comb_size(comb_matrix), 
                                                 ylim = c(0, max(comb_size(comb_matrix))*1.1),
                                                 add_numbers = TRUE,
                                                 border = FALSE,
                                                 gp = gpar(fill = "#5DA5DA", col = NA), 
                                                 height = unit(6, "cm")), 
            annotation_name_side = "left", 
            annotation_name_rot = 90,
            annotation_name_gp = gpar(fontsize = 8),
            annotation_name_offset = unit(1.2, "cm")),
          row_names_gp = gpar(fontsize = 8, fontface="bold"))

#Saving Upset plot
agg_png(file.path(xargs$outputDir,"UpSet_plot_for_genes_sets.png"), width = 4000, height = 3000, res = 600)
draw(f1,column_title = "Overlap of Differentially Expressed Genes (Ours vs. Authors)",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()


################# Comparison of logFC, baseMean, padj between gene sets of diffe-
################# rentially expressed genes
# "Up-regulated (Us only)",
# "Up-regulated (Authors only)",
# "Down-regulated (Us only)",
# "Down-regulated (Authors only)",
# "Up-regulated (Us & Authors)",
# "Down-regulated (Us & Authors)",
# "Up-regulated (Us) & Down-regulated (Authors)",
# "Down-regulated (Us) & Up-regulated (Authors)"

### Up-regulated genes from authors that aren't differentially expressed for us
# GenesID
ID.up.reg.genes.authors.not.DE.for.us = DE.genes.set$GeneID[DE.genes.set$up.regulated.authors==1 &
                                              DE.genes.set$down.regulated.authors==0 &
                                              DE.genes.set$up.regulated.us == 0 &
                                              DE.genes.set$down.regulated.us == 0]

# logFC, padj, baseMean from authors
up.reg.genes.authors.not.DE.for.us = differentially_expressed_genes.authors[differentially_expressed_genes.authors$GeneID
                                                                            %in% ID.up.reg.genes.authors.not.DE.for.us, 
                                                                            c("GeneID","baseMean","log2FoldChange","padj")]
up.reg.genes.authors.not.DE.for.us$From = rep("Authors", length(ID.up.reg.genes.authors.not.DE.for.us))
up.reg.genes.authors.not.DE.for.us$Status = rep("Up-regulated (Authors only)", length(ID.up.reg.genes.authors.not.DE.for.us))

# logFC, padj, baseMean from us
up.reg.genes.authors.not.DE.for.us.2 = dataf.DE.results.us[dataf.DE.results.us$GeneID %in% ID.up.reg.genes.authors.not.DE.for.us, 
                                                                            c("GeneID","baseMean","log2FoldChange","padj")]
up.reg.genes.authors.not.DE.for.us.2$From = rep("Us", length(ID.up.reg.genes.authors.not.DE.for.us))
up.reg.genes.authors.not.DE.for.us.2$Status = rep("Up-regulated (Authors only)", length(ID.up.reg.genes.authors.not.DE.for.us))

### Up-regulated genes from us that aren't differentially expressed for authors
# GenesID
ID.up.reg.genes.us.not.DE.for.authors = DE.genes.set$GeneID[DE.genes.set$up.regulated.authors==0 &
                                                              DE.genes.set$down.regulated.authors==0 &
                                                              DE.genes.set$up.regulated.us == 1 &
                                                              DE.genes.set$down.regulated.us == 0]

# logFC, padj, baseMean from us
up.reg.genes.us.not.DE.for.authors = differentially_expressed_genes.us[differentially_expressed_genes.us$GeneID
                                                                            %in% ID.up.reg.genes.us.not.DE.for.authors, 
                                                                            c("GeneID","baseMean","log2FoldChange","padj")]
up.reg.genes.us.not.DE.for.authors$From = rep("Us", length(ID.up.reg.genes.us.not.DE.for.authors))
up.reg.genes.us.not.DE.for.authors$Status = rep("Up-regulated (Us only)", length(ID.up.reg.genes.us.not.DE.for.authors))

# logFC, padj, baseMean from authors
up.reg.genes.us.not.DE.for.authors.2 = dataf.DE.results.authors[dataf.DE.results.authors$GeneID %in% ID.up.reg.genes.us.not.DE.for.authors, 
                                                           c("GeneID","baseMean","log2FoldChange","padj")]
up.reg.genes.us.not.DE.for.authors.2$From = rep("Authors", length(ID.up.reg.genes.us.not.DE.for.authors))
up.reg.genes.us.not.DE.for.authors.2$Status = rep("Up-regulated (Us only)", length(ID.up.reg.genes.us.not.DE.for.authors))


### Down-regulated genes from authors that aren't differentially expressed for us
# GenesID
ID.down.reg.genes.authors.not.DE.for.us = DE.genes.set$GeneID[DE.genes.set$up.regulated.authors==0 &
                                                              DE.genes.set$down.regulated.authors==1 &
                                                              DE.genes.set$up.regulated.us == 0 &
                                                              DE.genes.set$down.regulated.us == 0]

# logFC, padj, baseMean from authors
down.reg.genes.authors.not.DE.for.us = differentially_expressed_genes.authors[differentially_expressed_genes.authors$GeneID
                                                                            %in% ID.down.reg.genes.authors.not.DE.for.us, 
                                                                            c("GeneID","baseMean","log2FoldChange","padj")]
down.reg.genes.authors.not.DE.for.us$From = rep("Authors", length(ID.down.reg.genes.authors.not.DE.for.us))
down.reg.genes.authors.not.DE.for.us$Status = rep("Down-regulated (Authors only)", length(ID.down.reg.genes.authors.not.DE.for.us))

# logFC, padj, baseMean from us
down.reg.genes.authors.not.DE.for.us.2 = dataf.DE.results.us[dataf.DE.results.us$GeneID %in% ID.down.reg.genes.authors.not.DE.for.us, 
                                                           c("GeneID","baseMean","log2FoldChange","padj")]
down.reg.genes.authors.not.DE.for.us.2$From = rep("Us", length(ID.down.reg.genes.authors.not.DE.for.us))
down.reg.genes.authors.not.DE.for.us.2$Status = rep("Down-regulated (Authors only)", length(ID.down.reg.genes.authors.not.DE.for.us))


### Down-regulated genes from us that aren't differentially expressed for authors
# GenesID
ID.down.reg.genes.us.not.DE.for.authors = DE.genes.set$GeneID[DE.genes.set$up.regulated.authors==0 &
                                                              DE.genes.set$down.regulated.authors==0 &
                                                              DE.genes.set$up.regulated.us == 0 &
                                                              DE.genes.set$down.regulated.us == 1]

# logFC, padj, baseMean from us
down.reg.genes.us.not.DE.for.authors = differentially_expressed_genes.us[differentially_expressed_genes.us$GeneID
                                                                       %in% ID.down.reg.genes.us.not.DE.for.authors, 
                                                                       c("GeneID","baseMean","log2FoldChange","padj")]
down.reg.genes.us.not.DE.for.authors$From = rep("Us", length(ID.down.reg.genes.us.not.DE.for.authors))
down.reg.genes.us.not.DE.for.authors$Status = rep("Down-regulated (Us only)", length(ID.down.reg.genes.us.not.DE.for.authors))

# logFC, padj, baseMean from authors
down.reg.genes.us.not.DE.for.authors.2 = dataf.DE.results.authors[dataf.DE.results.authors$GeneID %in% ID.down.reg.genes.us.not.DE.for.authors, 
                                                                c("GeneID","baseMean","log2FoldChange","padj")]
down.reg.genes.us.not.DE.for.authors.2$From = rep("Authors", length(ID.down.reg.genes.us.not.DE.for.authors))
down.reg.genes.us.not.DE.for.authors.2$Status = rep("Down-regulated (Us only)", length(ID.down.reg.genes.us.not.DE.for.authors))


### Up-regulated genes from us that are down regulated for authors
# GenesID
ID.genes.up.us.down.authors = DE.genes.set$GeneID[DE.genes.set$up.regulated.authors==0 &
                                                                DE.genes.set$down.regulated.authors==1 &
                                                                DE.genes.set$up.regulated.us == 1 &
                                                                DE.genes.set$down.regulated.us == 0]

# logFC, padj, baseMean from us
genes.up.us.down.authors = differentially_expressed_genes.us[differentially_expressed_genes.us$GeneID
                                                                         %in% ID.genes.up.us.down.authors, 
                                                                         c("GeneID","baseMean","log2FoldChange","padj")]
genes.up.us.down.authors$From = rep("Us", length(ID.genes.up.us.down.authors))
genes.up.us.down.authors$Status = rep("Up-regulated (Us) & Down-regulated (Authors)", length(ID.genes.up.us.down.authors))

# logFC, padj, baseMean from authors
genes.up.us.down.authors.2 = differentially_expressed_genes.authors[differentially_expressed_genes.authors$GeneID %in% ID.genes.up.us.down.authors, 
                                                                  c("GeneID","baseMean","log2FoldChange","padj")]
genes.up.us.down.authors.2$From = rep("Authors", length(ID.genes.up.us.down.authors))
genes.up.us.down.authors.2$Status = rep("Up-regulated (Us) & Down-regulated (Authors)", length(ID.genes.up.us.down.authors))


### Down-regulated genes from us that are Up-regulated for authors
# GenesID
ID.genes.down.us.up.authors = DE.genes.set$GeneID[DE.genes.set$up.regulated.authors==1 &
                                                    DE.genes.set$down.regulated.authors==0 &
                                                    DE.genes.set$up.regulated.us == 0 &
                                                    DE.genes.set$down.regulated.us == 1]

# logFC, padj, baseMean from us
genes.down.us.up.authors = differentially_expressed_genes.us[differentially_expressed_genes.us$GeneID
                                                             %in% ID.genes.down.us.up.authors, 
                                                             c("GeneID","baseMean","log2FoldChange","padj")]
genes.down.us.up.authors$From = rep("Us", length(ID.genes.down.us.up.authors))
genes.down.us.up.authors$Status = rep("Down-regulated (Us) & Up-regulated (Authors)", length(ID.genes.down.us.up.authors))

# logFC, padj, baseMean from authors
genes.down.us.up.authors.2 = differentially_expressed_genes.authors[differentially_expressed_genes.authors$GeneID %in% ID.genes.down.us.up.authors, 
                                                                    c("GeneID","baseMean","log2FoldChange","padj")]
genes.down.us.up.authors.2$From = rep("Authors", length(ID.genes.up.us.down.authors))
genes.down.us.up.authors.2$Status = rep("Down-regulated (Us) & Up-regulated (Authors)", length(ID.genes.down.us.up.authors))

### Down-regulated genes from us that are Down-regulated for authors
# GenesID
ID.genes.down.us.down.authors = DE.genes.set$GeneID[DE.genes.set$up.regulated.authors==0 &
                                                    DE.genes.set$down.regulated.authors==1 &
                                                    DE.genes.set$up.regulated.us == 0 &
                                                    DE.genes.set$down.regulated.us == 1]

# logFC, padj, baseMean from us
genes.down.us.down.authors = differentially_expressed_genes.us[differentially_expressed_genes.us$GeneID
                                                             %in% ID.genes.down.us.down.authors, 
                                                             c("GeneID","baseMean","log2FoldChange","padj")]
genes.down.us.down.authors$From = rep("Us", length(ID.genes.down.us.down.authors))
genes.down.us.down.authors$Status = rep("Down-regulated (Us & Authors)", length(ID.genes.down.us.down.authors))

# logFC, padj, baseMean from authors
genes.down.us.down.authors.2 = differentially_expressed_genes.authors[differentially_expressed_genes.authors$GeneID %in% ID.genes.down.us.down.authors, 
                                                                    c("GeneID","baseMean","log2FoldChange","padj")]
genes.down.us.down.authors.2$From = rep("Authors", length(ID.genes.down.us.down.authors))
genes.down.us.down.authors.2$Status = rep("Down-regulated (Us & Authors)", length(ID.genes.down.us.down.authors))


### Up-regulated genes from us that are Up-regulated for authors
# GenesID
ID.genes.up.us.up.authors = DE.genes.set$GeneID[DE.genes.set$up.regulated.authors==1 &
                                                      DE.genes.set$down.regulated.authors==0 &
                                                      DE.genes.set$up.regulated.us == 1 &
                                                      DE.genes.set$down.regulated.us == 0]

# logFC, padj, baseMean from us
genes.up.us.up.authors = differentially_expressed_genes.us[differentially_expressed_genes.us$GeneID
                                                               %in% ID.genes.up.us.up.authors, 
                                                               c("GeneID","baseMean","log2FoldChange","padj")]
genes.up.us.up.authors$From = rep("Us", length(ID.genes.up.us.up.authors))
genes.up.us.up.authors$Status = rep("Up-regulated (Us & Authors)", length(ID.genes.up.us.up.authors))

# logFC, padj, baseMean from authors
genes.up.us.up.authors.2 = differentially_expressed_genes.authors[differentially_expressed_genes.authors$GeneID %in% ID.genes.up.us.up.authors, 
                                                                      c("GeneID","baseMean","log2FoldChange","padj")]
genes.up.us.up.authors.2$From = rep("Authors", length(ID.genes.up.us.up.authors))
genes.up.us.up.authors.2$Status = rep("Up-regulated (Us & Authors)", length(ID.genes.up.us.up.authors))



#### Full comparison table
DE.comparison.full.table = rbind(
  up.reg.genes.authors.not.DE.for.us,
  up.reg.genes.authors.not.DE.for.us.2,
  
  up.reg.genes.us.not.DE.for.authors,
  up.reg.genes.us.not.DE.for.authors.2,
  
  down.reg.genes.authors.not.DE.for.us,
  down.reg.genes.authors.not.DE.for.us.2,
  
  down.reg.genes.us.not.DE.for.authors,
  down.reg.genes.us.not.DE.for.authors.2,
  
  genes.up.us.down.authors,
  genes.up.us.down.authors.2,
  
  genes.down.us.up.authors,
  genes.down.us.up.authors.2,
  
  genes.down.us.down.authors,
  genes.down.us.down.authors.2,
  
  genes.up.us.up.authors,
  genes.up.us.up.authors.2
)

# Specifying the levels of Status
DE.comparison.full.table$Status <- factor(
  DE.comparison.full.table$Status,
  levels = c(
    "Up-regulated (Us only)",
    "Up-regulated (Authors only)",
    "Down-regulated (Us only)",
    "Down-regulated (Authors only)",
    "Up-regulated (Us & Authors)",
    "Down-regulated (Us & Authors)",
    "Up-regulated (Us) & Down-regulated (Authors)",
    "Down-regulated (Us) & Up-regulated (Authors)"
  )
)


############################ BOXPLOTS 

# For log2FC
f2 = ggplot(DE.comparison.full.table,
       aes(x = Status, y = log2FoldChange, fill = From)) +
  geom_boxplot(outlier.size = 0.5) +
  coord_flip() +
  labs(x = NULL,
       y = expression(log[2]~"Fold change"),
       fill = "From") +
  theme_bw()

ggsave(file.path(xargs$outputDir,"Boxplot_of_logFC_genes_set.png"), plot = f2, dpi = 600)


#For baseMean
f3 = ggplot(DE.comparison.full.table,
                 aes(x = Status, y = baseMean, fill = From)) +
  geom_boxplot(outlier.size = 0.5) +
  coord_flip() +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  labs(x = NULL,
       y = "Mean of normalized counts",
       fill = "From") +
  theme_bw()

ggsave(file.path(xargs$outputDir,"Boxplot_of_baseMeans_genes_set.png"), plot = f3, dpi = 600)

#For p-values
f4 = ggplot(DE.comparison.full.table,
            aes(x = Status, y = -log(padj, base = 10), fill = From)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  coord_flip() +
  labs(x = NULL,
       y = expression(-log[10]~"(padj)"),
       fill = "From") +
  theme_bw()

ggsave(file.path(xargs$outputDir,"Boxplot_of_adjusted_pvalue_genes_set.png"), plot = f4, dpi = 600)