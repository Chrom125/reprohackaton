#!/usr/bin/env Rscript
##################################### Plots ###############################


library(argparse) # For argument parsing
#Defining the arguments parser
parser = ArgumentParser(description= "This script performs a comparison between authors's differential expression analysis results and ours by generating boxplots and an UpSet plot.")


# Defining the arguments
parser$add_argument('--deseq2Results', '-dR', help= 'Path to DESeq2 results file', required= TRUE)
parser$add_argument('--outputDir', '-o', help= 'Path to the output directory', required= FALSE, default= './')


# Parsing the arguments
xargs = parser$parse_args()



###############################################################################
library(ggplot2)
library(scales) # For formatting axis graduations
library(ggrepel)# For better text label placement in ggplot2

##### MA plot for the entire RNAseq DataSet

dataf.DE.results = read.table(file = xargs$deseq2Results, header = T, sep = "\t", quote = "")

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
     scale_fill_manual(values = c("Non-Significant" = "black",
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

#===========> In KEGG database, S. aureus translation genes are grouped in 04 BRITE functional hierarchies
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
            aes(x = log(baseMean, base = 2) , y = log2FoldChange, color = signif)) +
  
     geom_point(alpha = 0.8, size = 0.9) + # Plotting all translation related genes
  
     # Plotting AA_tRNA_synthetases genes 
     geom_point(data = dataf.DE.results[filter_AA_tRNA_synthetases, ],
                shape = 1, color = "black", size = 1, stroke = 1.2) + 

     # Color scale for significant and non-significant genes
     scale_color_manual(values = c("Non-Significant" = "grey",
                                   "Significant" = "red")) +
     
     # Axis graduations
     scale_x_continuous(breaks = seq(0, 20, by = 2)) +
     scale_y_continuous(breaks = seq(-6, 5, by = 1)) +

     # Annotation segments for typical translation genes
     geom_segment(
      data = typical_members,
      aes(x = log(baseMean, base = 2) + 0.05, y = log2FoldChange,
          xend = log(baseMean, base = 2) + dx, yend = log2FoldChange + dy),
          linewidth = 0.8, color = "black") +
  
     # Annotation text for typical translation genes
     geom_text(
      data = typical_members,
      aes(x = log(baseMean, base = 2) + 1.35 * dx, y = log2FoldChange + 1.15 * dy, label = GeneName),
      inherit.aes = FALSE,
      fontface = "italic", size = 4) +
  
     # Horizontal line at y=0
     geom_hline(yintercept = 0, linetype = "dashed") +
  
     # Axis names + suppression of legend titles
     labs(
       x = expression(log[2]~"base Mean"),
       y = expression(log[2]~"fold change"), 
       color = NULL) +
     theme_classic() +
     theme(
      panel.border = element_rect(color = "black", linewidth = 2.5, fill = NA),
      
      # Putting the COLOR legend (Sig/Non-Sig) to the Bottom-Left (Inside)
      # Coordinates are relative (0 to 1)
      legend.position = c(0.12, 0.06), 
      legend.background = element_blank(), # Remove white box background
      legend.key = element_blank()         # Remove box around dots
     ) +
  
     # Manually creating the SHAPE legend (Synthetases) at Bottom-Right (Inside)
     annotate("point", x = 14, y = -5.5, shape = 1, size = 1, stroke = 1.2) +
     annotate("text", x = 14.5, y = -5.5, label = "AA-tRNA synthetases", 
              fontface = "italic", size = 3.5, hjust = 0)

    
# Saving
ggsave(file.path(xargs$outputDir,"MA_plot_translation_genes.png"), plot = f2, width = 6, height = 5, dpi = 300)