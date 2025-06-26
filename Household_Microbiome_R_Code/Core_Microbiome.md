# Core microbiome R code

Uses the .biom file from kraken-biom (kraken2 taxonomic analysis)

### Part 1 - setup

**Step 1 - Load required packages**

```
library(speedyseq)
library(phyloseq)
library(dplyr)
library(knitr)
library(microbiome)
library(microbiomeutilities)
library(stringr)
library(writexl)
```

**Step 2 - import biom file**

Import the biom file (in this example called biom_core) and fix metadata if needed (see Taxonomy_and_AlphaBetaDiv.md step 3.5 for instructions on this). 

```
biom_core <- import_biom(BIOMfilename = “PATH/TO/YOUR/FILE/biomfile.biom”)
```

### Part 2 - core microbiome calculations

**Step 1 - tranform .biom data**

Make sure data is **compositional** using the following command

```
biom_core.rel <- microbiome::transform(biom_core, "compositional")
```

**Step 2 - Aggregate taxa at different taxonomic levels**

You can change “Phylum” to any other level of interest

```
biom_phy_core <- aggregate_taxa(biom_core.rel, level="Phylum")
```

**Step 3 - calculate the core microbiome**

Find the core microbiome in the phyloseq object, where detection and prevelance can be changed to different numbers to specify how 'strict' the definition of core is (max=1, min=0)
Detection is the level at which each taxa must be found in each sample for it to be included (ie. detection = 0.0 means that if the taxa is included in a sample to ANY degree it will be included. Detection = 0.5 means that half of the sample must be that specific taxa in order to be included in the core). Prevalence is the number of samples which have the taxa in them (ie. prevalence = 0.99 means that 99% of the samples must have a taxa (at the level specified by detection) for that taxa to be included. Prevalance = 0.25 means that the taxa must be found in 25% of samples to be included). 

```
biom_core90.phy <- aggregate_rare(biom_phy_core, "unique", detection = 0.0, prevalence = 0.99)
biom_core90.phy <- taxa(biom_core90.phy)
print(lysol.all.core.taxa90.phy)
```

Data can be subsetted by different locations using the subset() function, and core microbiome analysis can be run on different sets of data (ie. all samples in one household or from one city) at different phylogenetic levels and with different detection and prevalence levels.

Note - the core microbiome for ARGs and VFs was determined using Venn diagrams (inputs were lists of the detected genes in each composite sample. Core genes were those found in the middle of the Venn diagram)

