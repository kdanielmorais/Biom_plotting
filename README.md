# Biom_plotting
Represent graphically the .biom data with .txt metadata

## Linux Tutorial

The program works as .R script, so in order to use it is necessary to activate it in the terminal.

The program will download all the library needed if necessary

> Rscript Biom_plotting.R 

The help function show the input arguments

> Rscript Biom_plotting.R --help

> -f CHARACTER, --file=CHARACTER
>		OTU table file name .BIOM [Required]

>	-t CHARACTER, --tax=CHARACTER
>		metadata file name .txt [Required]

>	-l CHARACTER, --level=CHARACTER
>		taxonomical level - choose between (Kingdom, Phylum, Class, Order, Family, Genus and Species) [default= Phylum]

>	-o CHARACTER, --out=CHARACTER
>		output file name [default= plot]

>	-h, --help
>		Show this help message and exit

After updating the input files the program will creat a .pdf file with the reespective ggplot graph

Into the files list there is a biom file and a metadata file for testing purposes




 
 

