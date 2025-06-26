# Global Household Microbiome

Example code for all analysis in the Global Household Microbiome paper. Code is written tutorial-style so it can easily be followed.

## ARG and VF analysis + taxonomic analysis
Code to identify antibiotic resistance genes, virulence factor genes, and to assign taxonomy to raw sequence data. Folder includes a README file with instructions/installation info, shell script to do the analysis, the specific code for each step in the analysis, example metadata for kraken2/generating a biom file, and two python scripts for sorting out ARG/VF genes from txt files/large excel sheets. This analysis was done on the University of Arizona's high performance computing system, but can be done on any computer with the resources to run the required software. 

## Household Microbiome R code
Code to run analysis VF/ARG/taxonomic data. Analysis was done using R in RStudio. Folder includes a README file with instructions/installation info, and .md files for each major analysis type (core microbiome, ARG and VF alpha and beta diversity, and taxonomic analysis with alpha and beta diversity). 

## Pathogen detection

Information on the python script (and the script itself) used to detect specific pathogens in the samples by screening kraken2 outputs.
