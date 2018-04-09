library("optparse")

cat ("                                                                                                                                                ")
cat ("Program version 1.0 --- GitHub link: https://github.com/ElcioNeto/Biom_plotting                                                                 ")
cat ("                                                                                                                                                ")
cat ("Exemples of program inputs:
     For a BIOM: Rscript SeedandBiom.R -c BIOM -m otu_table.biom -l Order -o plot_order
     For a SEED: Rscript SeedandBiom.R -c SEED -f otu_table.txt -t taxonomical_classification.txt -l Order -o plot_order                        ")
cat ("                                                                                                                                                ")
option_list = list(
  make_option(c("-c", "--cho"), type="character", dest= "choice", action= "store", default= "MAB",
              help="To utilize a biom file write BIOM and 
                to utilize a seed file write SEED [Required]"),
  make_option(c("-f", "--file"), type="character", default="MAV", 
              help="BIOM: a .biom Otu_table file
                SEED: a .txt Otu_table [Required]", metavar="character"),
  make_option(c("-t", "--tax"), type="character", default= "MAT", 
              help="SEED: a .txt taxonomical file [Required for seed]", metavar="character"),
  make_option(c("-m", "--met"), type="character", 
              help="A metadata file .txt [Optional for SEED and BIOM]", default= "MAC", metavar="character"),
  make_option(c("-l", "--level"),  type="character", default="Phylum", 
              help="taxonomical level - choose between (Kingdom, Phylum, Class, Order, Family, Genus and Species) [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="plot",
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

library("stringi")
library("biomformat")

## setting the with options
g <- opt$cho
h <- opt$met
a <- opt$file

if (g == "MAB") 
{
  cat("\n")
  cat("You have to choose a -c option,typing BIOM or SEED as entrance")
  cat("\n")
  
} else {

if (a == "MAV") 
  { 
  cat("\n")
  cat("You have to add a -f file for the reespective function chosen")
  cat("\n")
} else {
    
if (g == "BIOM") {
  
  ##  Getting the files added by the user and transforming it in a variablePegando os arquivos adicionados pelo usuário e transformando-os em variável

 
    ## dat <- read_biom(a)
  suppressWarnings(dat <- read_biom(a))
  
  otutable <- as.data.frame(as.matrix(biom_data(dat)))
  
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
  v <- 0
  for (v in 1:7) {
    
    taxonomy[,v] <- gsub("^$", "Unassigned", taxonomy[,v])
  }
 
  library("phyloseq")
  
  rltv = sweep(otucounts, 2, colSums(otucounts), "/")
  
  phy_otu = otu_table(rltv, taxa_are_rows = T)
  phy_taxa = tax_table(as.matrix(taxonomy))
  
  phy_total = phyloseq(phy_otu, phy_taxa)
  
  tax_level = opt$level
  
  cat("\n")
  cat(tax_level)
  cat("\n")
  
  phy_phy_glom = tax_glom(physeq = phy_total, taxrank = tax_level)
  
  df_melt = psmelt(physeq = phy_phy_glom)
  
  df_melt_aggregate_mean <- aggregate(df_melt$Abundance, by=list(Sample=df_melt$Sample, Taxa = df_melt[,tax_level]), FUN=sum)
  colnames(df_melt_aggregate_mean)<-c("Sample", tax_level, "Abundance")

  ## Plotting the graph 1
  
  library("ggplot2")
  
  #defining if plot legend separated or toguether
  tax_degree = nlevels(df_melt[,tax_level])
  
  ## In order to facilitate the visualization the legend and the graph will only be assembled in one single file if a small amount of taxonomical levels situation
  if(tax_degree <= 20) {
    pdf(paste("Sample_", opt$out,".pdf",sep = ""), width = 10, height = 10)
    plot = ggplot(data = df_melt_aggregate_mean, aes_string(x="Sample", y = "Abundance", fill=tax_level)) +
      geom_bar(stat = "identity", colour = "white")
    print(plot + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "relative abundance"))
    dev.off()
    print(paste("Samples plot was created in file ",getwd(),"/","Sample_mean", opt$out,".pdf", sep=""))
    
  }else 
  {
    pdf(paste("Sample_", opt$out,".pdf",sep = ""), width = 10, height = 10)
    plot = ggplot(data = df_melt_aggregate_mean, aes_string(x="Sample", y = "Abundance", fill=tax_level)) +
      geom_bar(stat = "identity", colour = "white")
    print(plot + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "relative abundance"))
    dev.off()
    
    
    #plotting legend isolated
    pdf(file = paste("legend_sample_",opt$out,".pdf",sep = ""), width = 12, height = 8)
    
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
    print("Two different .pdf files of samples were created separatedly because of the high number of taxonomical levels")
  }
  
  cat("\n")
  
  ## Default samples graph
  if (h == "MAC") 
  {
    cat("\n")
    cat("Without a metadata file is only possible to create a plot of the relative abundance x samples ")
    cat("\n")
  } else {
  
  ## Processing with a group mean utilizing a metadata file  
  
  ## obtaining the metadata file  
  met <- read.delim(h)
  
  ## Adapting the format of the metadata to be utilized
  rownames(met) = met[,1]
  
  ## Using the phyloseq function to the data grouping
  phy_data = sample_data(met)
  
  phy_total1 = phyloseq(phy_otu, phy_data, phy_taxa) 
  
  ## Grouping according to the choiced taxonomical level
  phy_phy_glom1 = tax_glom(physeq = phy_total1, taxrank = tax_level)
  
  ## melting the file for porterior use of ggplot function 
  df_melt1 = psmelt(physeq = phy_phy_glom1)
  
  ## Melting the file for the group mean analysis
  df_melt_aggreg = aggregate(df_melt1$Abundance, by=list(Sample = df_melt1$Sample, Group=df_melt1$Group, taxa=df_melt1[,tax_level]), FUN=sum)
  colnames(df_melt_aggreg)<-c("Sample","Group", tax_level, "Abundance")
  
  df_melt_aggreg_mean1 = aggregate(df_melt_aggreg$Abundance, by=list(Group=df_melt_aggreg$Group, taxa=df_melt_aggreg[,tax_level]), FUN=mean)
  colnames(df_melt_aggreg_mean1)<-c("Group", tax_level, "Abundance")
  
  ## Plotting the graph 2
  
  ## In order to facilitate the visualization the legend and the graph will only be assembled in one single file if a small amount of taxonomical levels situation
  if(tax_degree <= 20) {
    pdf(paste("Group_", opt$out,".pdf",sep = ""), width = 10, height = 10)
    plot = ggplot(data = df_melt_aggreg_mean1, aes_string(x="Group", y = "Abundance", fill=tax_level)) +
      geom_bar(stat = "identity", colour = "white")
    print(plot + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "relative abundance"))
    dev.off()
    print(paste("Group plot was created in file ",getwd(),"/", "Group_mean_", opt$out,".pdf", sep=""))
    
  }else 
  {
    pdf(paste("Group_", opt$out,".pdf",sep = ""), width = 10, height = 10)
    plot = ggplot(data = df_melt_aggreg_mean1, aes_string(x="Group", y = "Abundance", fill=tax_level)) +
      geom_bar(stat = "identity", colour = "white")
    print(plot + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "relative abundance"))
    dev.off()
   
     #plotting legend isolated
    pdf(file = paste("legend_group_",opt$out,".pdf",sep = ""), width = 12, height = 8)
    
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
    print("Two different .pdf files of relative abundance-x-group  were created separatedly because of the high number of taxonomical levels")
  }
  }
}

if (g == "SEED") 
{
  
  d <- opt$tax
 
  if (d == "MAT") 
    {
     cat("Its not possible to process the SEED graph without a taxonomical file")
     cat("\n")
  } else{
  
  df_otu_table = read.delim(a)
  tax_table = read.delim(d)
  
  ##Parsing data
  #Table of counts
  data = df_otu_table[(-nrow(df_otu_table)),2:(ncol(df_otu_table)-2)]
  rownames(data) = data[,1]
  otucounts = data[,-1]
  
  #Tax table
  df = tax_table[,4]
  x = strsplit(as.character(df), "\\|")
  x2 = lapply(x, FUN=function(df) unlist(strsplit(df,";")))
  taxonomy_info = as.data.frame(stri_list2matrix(x2, byrow = T, fill = "unassigned"))
  colnames(taxonomy_info) = c("Unite_SpHip","NCBI_accession","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  rownames(taxonomy_info) = data[,1]
  
  library("phyloseq")
  
  #generating relative abundance counts
  rltv = sweep(otucounts, 2, colSums(otucounts), "/")
  
  #generating a phyloseq objects
  rtlOTU = otu_table(rltv, taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(taxonomy_info))
  rltfungiPhy = phyloseq(rtlOTU, TAX)
  
  #melting phyloseq objects for ggplot manipulation
  rltmds = psmelt(rltfungiPhy)
  
  tax_level = opt$level
  cat("\n")
  cat(tax_level)
  cat("\n")
  
  rltmds_aggregate <- aggregate(rltmds$Abundance, by=list(Sample=rltmds$Sample, Taxa = rltmds[,tax_level]), FUN=sum)
  colnames(rltmds_aggregate)<-c("Sample", tax_level, "Abundance")
  
  ##Plotting
  library("ggplot2")
  
  #defining if plot legend separated or toguether
  tax_degree = nlevels(rltmds[,tax_level])
  
  if(tax_degree <= 20) {
    pdf(paste("Sample_", opt$out,".pdf",sep = ""), width = 10, height = 10)
    plot = ggplot(data = rltmds_aggregate, aes_string(x="Sample", y = "Abundance", fill = tax_level)) +
      geom_bar(stat = "identity", colour = "white")
    print(plot + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "relative abundance"))
    dev.off()
    print(paste("Sample plot was created in file ",getwd(),"/", "Group_", opt$out,".pdf", sep=""))
  } else {
    pdf(paste("Sample_mean_", opt$out,".pdf",sep = ""), width = 10, height = 10)
    plot = ggplot(data = rltmds_aggregate, aes_string(x="Sample", y = "Abundance", fill = tax_level)) +
      geom_bar(stat = "identity", colour = "white")
    print(plot + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "relative abundance"))
    dev.off()
    
    #plotting legend isolated
    pdf(file = paste("legend_sample_",opt$out,".pdf",sep = ""), width = 12, height = 8)
    
    plot_legend <- function(a.gplot){
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
    print("Sample legend and Plot were created in separate files in reason of the high number of taxonomical levels")
  }
  
  if (h == "MAC") 
    {
    cat("\n")
    cat("Without a metadata file is only possible to create a plot of the relative abundance x samples")
    cat("\n")
  }else 
    {
    ## Processing with a group mean utilizing a metadata file  
      
    ## obtaining the metadata file  
    met <- read.delim(b)
      
    ## Using the phyloseq function to the data grouping
    phy_data = sample_data(met)
    rltfungiPhy1 = phyloseq(rtlOTU, phy_data, TAX)
    
    #melting phyloseq objects for ggplot manipulation
    rltmds1 = psmelt(rltfungiPhy1)
    
    ## Melting the file for the group mean analysis
    rltmds_aggregate1 <- aggregate(rltmds1$Abundance, by=list(Sample=rltmds1$Sample, Group=rltmds1$Group, Taxa = rltmds1[,tax_level]), FUN=sum)
    colnames(rltmds_aggregate)<-c("Sample", "Group", tax_level, "Abundance")
    
    rltmds_aggregate2 <- aggregate(rltmds_aggregate1$Abundance, by=list(Group=rltmds_aggregate1$Group, Taxa = rltmds_aggregate1[,tax_level]), FUN=mean)
    colnames(rltmds_aggregate2)<-c("Group", tax_level, "Abundance")
    
    ## plotting graph 2
    
    #defining if plot legend separated or toguether
    tax_degree = nlevels(rltmds[,tax_level])
    
    if(tax_degree <= 20) {
      pdf(paste("Group_", opt$out,".pdf",sep = ""), width = 10, height = 10)
      plot = ggplot(data = rltmds_aggregate2, aes_string(x="Group", y = "Abundance", fill = tax_level)) +
        geom_bar(stat = "identity", colour = "white")
      print(plot + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "relative abundance"))
      dev.off()
      print(paste("Group plot was created in file ",getwd(),"/", "Group_", opt$out,".pdf", sep=""))
    } else {
      pdf(paste("Group_", opt$out,".pdf",sep = ""), width = 10, height = 10)
      plot = ggplot(data = rltmds_aggregate2, aes_string(x="Group", y = "Abundance", fill = tax_level)) +
        geom_bar(stat = "identity", colour = "white")
      print(plot + theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "relative abundance"))
      dev.off()
      
      #plotting legend isolated
      pdf(file = paste("legend_group_",opt$out,".pdf",sep = ""), width = 12, height = 8)
      
      plot_legend <- function(a.gplot){
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
      }
      #Plotting only legendss. Important for species and genus factors. Too many levels.
      mylegend = plot_legend(plot)
      library(grid)
      grid.draw(mylegend)
      dev.off()
      print("Two .pdf files of relative abundance-x-group were created in separate files because of the high number of taxonomical levels")
      
    }  
    }  
    }
    }
    }
  }
  
  
  
 