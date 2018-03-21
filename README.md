# Biom and Seed plotting
Represent graphically the relative abundance of a .biom or .txt data using the samples as the second parameter or the groups if the user add a metadata file.

## Linux Tutorial

The program works as .R script, so in order to use it is necessary to activate it in the terminal.

> Rscript BiomandSeed_plotting_end.R 

The help function show the input arguments

> Rscript Biom_plotting.R --help

> -c CHO, --cho=CHO
		To utilize a biom file write BIOM and 
  to utilize a seed file write SEED [Required]

	-f CHARACTER, --file=CHARACTER
		BIOM: a .biom Otu_table file
                SEED: a .txt Otu_table [Required]

	-t CHARACTER, --tax=CHARACTER
		SEED: a .txt taxonomical file [Required for seed]

	-m CHARACTER, --met=CHARACTER
		A metadata file .txt [Optional for SEED and BIOM]

	-l CHARACTER, --level=CHARACTER
		taxonomical level - choose between (Kingdom, Phylum, Class, Order, Family, Genus and Species) [default= Phylum]

	-o CHARACTER, --out=CHARACTER
		output file name [default= plot]

	-h, --help
		Show this help message and exit


After updating the input files the program will creat a .pdf file with the reespective ggplot graph

Into the files list there is a biom file and a metadata file for testing purposes




 
 

