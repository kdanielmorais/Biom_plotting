# Biom and Seed plotting
The script represents graphically the relative abundance of a .biom or .txt data using the samples as the second parameter or the groups if the user add a metadata file.
There is one Rscript for the 16S database and one for the ITS database

## Tutorial

The program works as .R script, so in order to use it is necessary to access his diretory and activate it in the terminal.

> Rscript BiomandSeed_plotting_16S.R 

The help function show the  following input arguments

> Rscript BiomandSeed_plotting_16S.R --help

 	-c CHO, --cho=CHO
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


Its necessary to add at least -c and -f entrances in order to the script produce the graph. If those entrances do not be filled the script will ask for them.

If the user do not choose the level, the script chooses the option "Phylum".

If the user do not choose the output name, the script chooses the name "plot", so it will generate a Sample_plot file and if the option metadata be filled it will create a Group_plot file as well.

In the -f input file has a great amount of taxonomical levels the script will generate two differente .pdf files for each graph, one with the graph itself and another with the legend. The legend file will have the same name as the graph but with "Legend" inserted. 

The BIOM analysis requires only the .biom file to generate the sample graph

The SEED analysis requires one otu table file .txt and one taxonomical file .txt. If the taxonimical file do not be fullfiled the script will ask for it

Each time the script is used he reproduces the version and his github link. He will reproduce two following entrance exemples as well:

> For a BIOM: Rscript BiomandSeed_plotting_16S.R -c BIOM -m otu_table.biom -l Order -o plot_order
> For a SEED: Rscript BiomandSeed_plotting_16S.R -c SEED -f otu_table.txt -t taxonomical_classification.txt -l Order -o plot_order 

In the archives setted in Github there is a biom file, a metadata for the biom, a otu table .txt for seed and a taxonomic for seed to be used as examples.

There is the examples of samples and group graphs to show the appearence of the graphs generated.
