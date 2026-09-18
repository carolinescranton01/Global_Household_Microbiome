# Core microbiome R code

Uses the decontaminated biomfile object from the taxonomic and alpha/beta diversity analysis (https://github.com/carolinescranton01/cs_projects/blob/main/Household_Microbiome_R_Code/Taxonomy_and_AlphaBetaDiv.md) 


### Part 1 - setup

**Step 1 - Load required packages**

```
library(phyloseq) # Version 1.52.0
library(dplyr) # Version 1.2.0
library(knitr) # Version 1.51
library(microbiome) # Version 1.30.0
library(microbiomeutilities) # Version 1.0.17
library(stringr) # Version 1.6.0
library(writexl) # Version 1.5.4
```

**Step 2 - import biom file**

Import the biomfile object (in this example called biom_core) and fix metadata if needed (see Taxonomy_and_AlphaBetaDiv.md step 3.5 for instructions on this). 

```
biom_core <- biomfile
print(biom_core)
```

# Create a new column for full organism names (Genus + Species), update phyloseq object, transform to relative abundance

```
tax_table_full <- as.data.frame(tax_table(biom_core))
tax_table_full$FullName <- paste(tax_table_full$"Genus", tax_table_full$"Species")

# update the tax_table in the phyloseq object

tax_table(biom_core) <- as.matrix(tax_table_full)

core_rel <- microbiome::transform(biom_core, "compositional")
```

### Part 2 - core microbiome calculations

**Step 1 - aggregate at taxonomic levels**

Make sure data is **compositional**, and aggregate at different taxonomic levels for different analyses. You can change “Phylum” to any other level of interest

```
core_rel.phy <- aggregate_taxa(core_rel, level="Phylum")
```

**Step 2 - calculate the core microbiome**

Find the core microbiome in the phyloseq object, where detection and prevelance can be changed to different numbers to specify how 'strict' the definition of core is (max=1, min=0)
Detection is the level at which each taxa must be found in each sample for it to be included (ie. detection = 0.0 means that if the taxa is included in a sample to ANY degree it will be included. Detection = 0.5 means that half of the sample must be that specific taxa in order to be included in the core). Prevalence is the number of samples which have the taxa in them (ie. prevalence = 0.99 means that 99% of the samples must have a taxa (at the level specified by detection) for that taxa to be included. Prevalance = 0.25 means that the taxa must be found in 25% of samples to be included). 

```
biom_core90.phy <- aggregate_rare(biom_phy_core, "unique", detection = 0.0, prevalence = 0.99)
biom_core90.phy <- taxa(biom_core90.phy)
print(lysol.all.core.taxa90.phy)
```

Data can be subsetted by different locations using the subset() function, and core microbiome analysis can be run on different sets of data (ie. all samples in one household or from one city) at different phylogenetic levels and with different detection and prevalence levels.

Note - the core microbiome for ARGs and VFs was determined using Venn diagrams (inputs were lists of the detected genes in each composite sample. Core genes were those found in the middle of the Venn diagram)

