
##Checking packages (CRAN)
list.of.packages = c("optparse", "ggplot2", "stringi", "BIOM.utils")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages())]
if(length(new.packages)) 
{
  install.packages(new.packages,  repos = "http://cran.us.r-project.org")
}

##Checking Packages (Bioconductor)
list.of.packages = c("phyloseq", "biomformat")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages())]
if(length(new.packages)) 
{
  source('http://bioconductor.org/biocLite.R')
  biocLite('phyloseq' , "biomformat")
 
}

## Setting the library
library("stringi")
library("phyloseq")
library("ggplot2")
library("BIOM.utils")
library("biomformat")
library("optparse")


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="OTU table file name .BIOM [Required]", metavar="character"),
  make_option(c("-t", "--tax"), type="character", default=NULL,
              help="metadata file name .txt [Required]", metavar="character"),
  make_option(c("-l", "--level"),  type="character", default="Phylum", 
              help="taxonomical level - choose between (Kingdom, Phylum, Class, Order, Family, Genus and Species) [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="plot",
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least two input files must be supplied. \n--file [OTU_table.txt] and --tax [metadata_table.txt] )\n", call.=FALSE)
}

##  Getting the files added by the user and transforming it in a variablePegando os arquivos adicionados pelo usuário e transformando-os em variável
a <- opt$file
b <- opt$tax

dat <- read_biom(a)
met <- read.table(b)

## Transforming the .biom file in an OTU table
otutable <- as.data.frame(as.matrix(biom_data(dat)))

## Adapting the format of the metadata to be utilized
rownames(met) = met[,1]
colnames(met) = c("Sample ID","Barcode Sequence","Linker Primer Sequence","Group","Description")

## Getting the taxonomy
taxonomy3 <- observation_metadata(dat)

## Renaming the OTU
df_otu_table <- otutable
otucounts <- otutable

## Removing the column "confidence" from the taxanomy 
taxonomy <- taxonomy3[,-c(1)]
tax_table <- taxonomy

## Organizing the table and changing the titles of the columns
colnames(taxonomy) = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
taxonomy_info1 = as.data.frame(taxonomy, byrow = T, fill = "unassigned")
taxonomy <- taxonomy_info1


## Generating the relative abundance
rltv = sweep(otucounts, 2, colSums(otucounts), "/")

## Using the phyloseq function to the data grouping
phy_otu = otu_table(rltv, taxa_are_rows = T)
phy_data = sample_data(met)
phy_taxa = tax_table(as.matrix(taxonomy))

phy_total = phyloseq(phy_otu, phy_data, phy_taxa) 

## Selecting the taxonomical level 
tax_level = opt$level

## Grouping according to the choiced taxonomical level
phy_phy_glom = tax_glom(physeq = phy_total, taxrank = tax_level)
colSums(otu_table(phy_phy_glom)) 

## melting the file for porterior use of ggplot function 
df_melt = psmelt(physeq = phy_phy_glom)
colnames(df_melt)

df_melt_aggreg1 = aggregate(df_melt$Abundance, by=list(Sample = df_melt$Sample, Group=df_melt$Group, taxa=df_melt[,tax_level]), FUN=sum)
colnames(df_melt_aggreg1)<-c("Sample","Group", tax_level, "Abundance")

df_melt_aggreg_mean = aggregate(df_melt_aggreg1$Abundance, by=list(Group=df_melt_aggreg1$Group, taxa=df_melt_aggreg1$Phylum), FUN=mean)
colnames(df_melt_aggreg_mean)<-c("Group", tax_level, "Abundance")

## Plotting the graph

#defining if plot legend separated or toguether
tax_degree = nlevels(df_melt[,tax_level])

## In order to facilitate the visualization the legend and the graph will only be assembled in one single file if a small amount of taxonomical levels situation
if(tax_degree <= 20) {
  pdf(paste(opt$out,".pdf",sep = ""), width = 10, height = 10)
  plot = ggplot(data = df_melt_aggreg_mean, aes_string(x="Group", y = "Abundance", fill=tax_level)) +
    geom_bar(stat = "identity", colour = "white")
  print(plot + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "relative abundance"))
  dev.off()
  print(paste("Plot was created in file ",getwd(),"/",opt$out,".pdf", sep=""))
  
}else 
{
  pdf(paste(opt$out,".pdf",sep = ""), width = 10, height = 10)
  plot = ggplot(data = df_melt_aggreg_mean, aes_string(x="Group", y = "Abundance", fill=tax_level)) +
    geom_bar(stat = "identity", colour = "white")
  print(plot + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "relative abundance"))
  dev.off()
  
  
  #plotting legend isolated
  pdf(file = paste("legend_",opt$out,".pdf",sep = ""), width = 12, height = 8)
  
  plot_legend <- function(a.gplot)
  {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  #Plotting only legend. Important for species and genus factors. Too many levels.
  mylegend = plot_legend(plot)
  library(grid)
  grid.draw(mylegend)
  dev.off()
  print("Two different .pdf files were created. Legend and Plot were created in separate files in reason of the high number of taxonomical levels")
}
