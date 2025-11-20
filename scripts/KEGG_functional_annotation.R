############################# KEGG functional annotation ################################

library(argparse) # For parsing command line arguments

#Defining the arguments parser
parser = ArgumentParser(description= 'This scripts downloads the KEGG BRITE functional annotation for S.aureus NCTC 8325 genes 
and saves it in a tsv file.')

# Defining the arguments
parser$add_argument('--outputDir', '-o', help= 'Path to the output directory', 
                    required= FALSE, default= './')

# Parsing the arguments
xargs = parser$parse_args()


library(KEGGREST) # For accessing KEGG database : KEGG REST API

#### Getting KEGG BRITE functional hierarchy info specifying metabolic pathways for S.aureus (strain NCTC 8325) genes
#### S. aureus code in KEGG database is "sao"

#Correspondance table between S. aureus (code sao) genes and their KEGG BRITE functional hierarchy ID
GeneID_BRITE.ID= keggLink(target = "brite", source = "sao") # vector of KEGG BRITE functional hierarchy (with the genes as names)

GeneID = sub(pattern = "^sao:",  replacement = "", names(GeneID_BRITE.ID)) #genes ID without the prefix "sao:"

BRITE.ID = sub(pattern = "^br:",  replacement = "", GeneID_BRITE.ID) #KEGG BRITE functional hierarchy ID without the prefix "br:"

GeneID_BRITE.ID= data.frame(GeneID, BRITE.ID) # table GeneID | KEGG BRITE ID


#Aggregating all the hierarchies for each gene, since a gene may be active in many functional pathway
GeneID_BRITE.ID = aggregate(. ~ GeneID, data = GeneID_BRITE.ID,
                             FUN = function(x) paste(unique(x), collapse = ";"))

#Writing the functional annotation table in a tsv file
write.table(GeneID_BRITE.ID, file = file.path(xargs$outputDir,"KEGG_BRITE_functional_annotation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)