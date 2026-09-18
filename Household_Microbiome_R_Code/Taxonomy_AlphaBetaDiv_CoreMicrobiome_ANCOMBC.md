## RMarkdown walkthrough on how to analyze .biom files for taxonomic composition, alpha diversity, and beta diversity

This code was used to analyze all samples for both bacterial, viral, and eukaryotic pathogen taxa, with a few minor edits (such as changing the taxonomic ranks for viruses)

### Part 1 - setup and data cleanup

**Step 1. Load required packages**

```
# note - you may need a package manager. BiocManager works well. 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("microbiome") # Version 1.30.0
BiocManager::install("phyloseq") # Version 1.52.0
BiocManager::install("microbiomeutilities") # Version 1.0.17
BiocManager::install("RColorBrewer") # Version 1.1.3
BiocManager::install("ggpubr") # Version 0.6.3
BiocManager::install("DT") # Version 0.34.0
BiocManager::install("data.table") # Version 1.18.2.1
BiocManager::install("dplyr") # Version 1.2.0
BiocManager::install("writexl") # Version 1.5.4
BiocManager::install("openxlsx") # Version 4.8.2.1
BiocManager::install("vegan") # Version 2.7.3
BiocManager::install("ggplot2") # Version 4.0.2
BiocManager::install("decontam") # Version 1.32.0
BiocManager::install("multcompView") # Version 0.1.11
BiocManager::install("stringr") # Version 1.6.0
BiocManager::install("ANCOMBC") # Version 2.14.0


# If packages are already installed, load via library()

library(microbiome) 
library(phyloseq) 
library(microbiomeutilities) 
library(RColorBrewer)
library(ggpubr)
library(DT)
library(data.table)
library(dplyr)
library(writexl)
library(openxlsx)
library(vegan)
library(ggplot2)
library(decontam)
library(multcompView)
library(stringr)
library(ANCOMBC)
```

**Step 2: Import biom file:**
```
biomfile <- import_biom(BIOMfilename = “PATH/TO/YOUR/FILE/biomfile.biom”)
```

**Step 3: Change column names in taxonomy to taxonomic levels (instead of ‘Rank 1’, ‘Rank 2’, etc)**

For bacteria, use c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"). For viruses, use c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

```
colnames(tax_table(biomfile)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
```

**Step 3.5 (OPTIONAL): If needed, replace metadata within biom file using excel file with updated metadata (has to have same Sample_ID as biom file)**

If you need to add additional metadata to the .biom file, create a new excel file with the entire (updated) metadata set. Make sure the sample_ID column matches the sample_ID column in the original metadata/the .biom file - the actual values within the column as well as the header must match for this to work! In this example, the column header is Sample_ID.

```
updated_metadata_1 <- read.xlsx("updated_metadata.xlsx")
rownames(updated_metadata_1) <- updated_metadata_1$Sample_ID
updated_metadata_2 <- sample_data(updated_metadata_1)
updated_metadata_2$Sample_ID <- NULL
biomfile <- merge_phyloseq(biomfile, updated_metadata_2)
```

**Step 4: Remove unwanted taxa**

The database used to assign taxonomic IDs for this dataset included all known and non-redundant sequences in the NCBI database (for the viral analysis, this is not the case - the database used only included viral DNA). We need to remove some of these sequences as they are considered contaminants in the dataset.

```
#remove mitochondira
biomfile_nomitochondria <- subset_taxa(biomfile, Family != "Mitochondria")
ntaxa(biomfile)-ntaxa(biomfile_nomitochondria)
#remove chloroplast
biomfile_nochloroplast<- subset_taxa(biomfile, Family != "Chloroplast")
ntaxa(biomfile)-ntaxa(biomfile_nochloroplast)
# remove Human DNA
biomfile_nohuman <- subset_taxa(biomfile, Family != "Hominidae")
ntaxa(biomfile)-ntaxa(biomfile_nohuman)
biomfile_nohuman2 <- subset_taxa(biomfile, Genus != "Homo")
ntaxa(biomfile)-ntaxa(biomfile_nohuman2)

# Alternativly, filter to only include kingdom Bacteria
biomfile <- subset_taxa(biomfile, Domain == "Bacteria")
```

Now we should just be left with bacteria reads

Decontaminate reads using decontam package and negative control samples
```
# create a new columns for negative controls, based on metadata value. In this case there was a column called Sample_Type, with values "Negative_Control" and "Sample". This command makes a new column for this analysis called is.neg, where the values are either TRUE or FALSE
sample_data(biomfile)$is.neg <-
  sample_data(biomfile)$Sample_Type == "Negative_Control"

# this identifies contaminants based on the prevalence of taxa in the negative control samples vs the test samples
contam_df <- isContaminant(
  biomfile,
  method = "prevalence",
  conc = "reads",
  neg="is.neg"
)

# this creates a table which lists the contaminant taxa
contaminants_to_remove <- rownames(contam_df)[contam_df$p <= 0.05]

# this removes contaminants from the sample data
biomfile_clean <- prune_taxa(
  !(taxa_names(biomfile) %in% contaminants_to_remove),
  biomfile
)

# generate read statistics before/after contaminant filtering as a quality check

# total number of reads per sample before decontamination
raw_reads <- sample_sums(biomfile_clean)
# total number of reads per sample after decontamination
clean_reads <- sample_sums(biomfile)

# summarize the total reads and reads lost
read_summary <- tibble(
  Status = c("Before Decontam", "After Decontam"),
  Total_Reads = c(sum(raw_reads), sum(clean_reads)),
  Average_Reads = c(mean(raw_reads), mean(clean_reads)),
  Min_Reads = c(min(raw_reads), min(clean_reads)),
  Max_Reads = c(max(raw_reads), max(clean_reads))
)
print(read_summary)

# count the number of taxa lost
ntaxa_before <- ntaxa(biomfile)
ntaxa_after <- ntaxa(biomfile_clean)
taxa_lost <- ntaxa_before - ntaxa_after

message(paste("Taxa before decontam:", ntaxa_before))
message(paste("Taxa after decontam:", ntaxa_after))
message(paste("Total Taxa removed:", taxa_lost))

# OPTIONAL - rename biomfile_clean to biomfile for consistence with downstream code, remove negative control samples (no longer needed)
biomfile <- biomfile_clean

# remove negative control samples
biomfile  <- subset_samples(
  biomfile,
  Sample_Type != "Negative_Control"
)
```


### Part 2 - figures and diversity analysis

**Preliminary figures to look at phyla prevalence, sequencing depth, etc**

These figures are generated to assess the data quality and completeness overall. 

```
#Phylum plot
biom_1 <- plot_taxa_cv(biomfile, plot.type = "scatter")
biom_1 + scale_x_log10()

#Sequencing depth by geographic location - can be changed to any column in the metadata
biom_seqdepth.ngtax_Geo <- plot_read_distribution(biomfile, "Geographic_Location", "density")
print(biom_seqdepth.ngtax_Geo)

#Histogram of ASVs – reformat data structure and then graph
biom_histogram_data<- data.table(
  	tax_table = as.data.frame(tax_table(biomfile)),
  	ASVabundance = taxa_sums(biomfile),
  	ASV = taxa_names(biomfile))
biom_histogram_plot <- ggplot(biom_histogram_data, aes(ASVabundance)) +  
geom_histogram() +
    ggtitle("Histogram of ASVs (unique sequence) counts") +
    theme_bw() +
    scale_x_log10() +
    ylab("Frequency of ASVs") +
    xlab("Abundance (raw counts)")
print(biom_histogram_plot)
```

**Rarefication of samples + plot for beta diversity analysis**

Samples must be rarefied to account for some sequences having less data than others

```
set.seed(1234)
biom_rar <- rarefy_even_depth(biomfile, sample.size = 10000) # choose as high as possible number that does not lead to significant data loss
print(biom_rar)
barplot(sample_sums(biom_rar), las =2)
```

**Alpha diversity calculations**

```
biom.alphadiv <- alpha(biom_rar, index = "all")
biom.alphadiv$SampleID <- rownames(biom.alphadiv) # sample IDs are the last column
diversity_dataframe <- merge(biom.alphadiv, metadata, by = "SampleID") # note: change metadata to your metadata object name if needed
```

Below is an example of the code for a plot with statistics. X can be changed to different columns in the metadata, and Y can be changed to different alpha diversity metrics, calculated with the command above. Tukey's post-hoc test can be done for significant comparisons

```
household_shannon_div <- ggboxplot(diversity_dataframe, x = "Household_Location", y = "diversity_shannon") +
    rotate_x_text() +
    ylim(0, 10) +
    theme(legend.position="none") +
    labs(title="Shannon Diversity in Households in Different Cities") +
    stat_compare_means(aes(group=Household_Location), label = 'p.format', method='anova', label.y=5, label.x=1.5) +
    facet_wrap(~Geographic_Location)

# Tukey's post-hoc test
house_aov <- aov(diversity_shannon ~ Household_Location, data = diversity_dataframe)
tuk <- TukeyHSD(house_aov)
letters <- multcompLetters4(house_aov, tuk)
print(letters) # prints significance letters, where variables that share letters (ie a, ab, and abc) are not significantly different from each other
```

**Beta diversity calculations and PERMANOVA test**

PCoA ordination, PERMANOVA, and plot from phyloseq object. Data can be subsetted before ordination to run an analysis on specific sets of samples

```
dist_bray <- phyloseq::distance(biom_rar, method = "bray")
metadata <- data.frame(sample_data(biom_rar))
PERMANOVA_result <- adonis2(dist_bray ~ Geographic_Location, data = metadata)
print(PERMANOVA_result)

biom_data <- ordinate(biom_rar, "PCoA", "bray")
plot_ordination(biom_rar, biom_data, "bray", color = "Geographic_Location", shape = "Geographic_Location") +
    geom_point(size = 1) +
    ggtitle("PCoA of Data") +
    font("ylab", size = 12, face = "bold") + 
    stat_ellipse(aes(color = Geographic_Location), level = 0.95, size = 0.5) + 
    font("xlab", size = 12, face = "bold") + font("title", size = 10, face = "bold")
```

**Generating taxonomic relative abundance graphs**

```
biom_phylum <- aggregate_top_taxa2(biomfile, top = 10, "Phylum") 
biom_phylum.rel <- microbiome::transform(biom_phylum, "compositional")
dfphy <- psmelt(biom_phylum.rel)

# Convert phyloseq object to a data frame
biom_phy.rel.abun <- ggplot(dfphy, aes(x = Household_Location, y = Abundance, 
fill = Phylum)) +
geom_bar(stat = "identity") +
 	facet_wrap(~ Geographic_Location, scales = "free_x") +
  theme(legend.position = "bottom", legend.text = element_text(size = 6), legend.title = element_text(size = 10), legend.key.size = unit(0.5, "cm"), 
  legend.spacing.y = unit(0.1, "cm")) +
 	scale_fill_brewer("Phylum", palette = "Paired") + 
 	theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, size = 6)) + 
  labs(title = "Relative Abundance by Household Location", x = "Household Location", y = "Relative Abundance") +
  guides(fill = guide_legend(title = "Phylum", title.theme = element_text(size = 10), label.theme = element_text(size = 8), keywidth = unit(0.5, "cm"), keyheight = unit(0.5, "cm")))

print(biom_phy.rel.abun)
```

**NOTE:** The same code as above were used for other taxonomic levels. 'X' was changed to geographic_location when looking at each city overall. Data can be subsetted by taxonomic ranks to look at specific geographic or household location’s taxonomy

### Part 3 - core microbiome calculations

**Step 1 - import biom file**

Use the same biomfile object as previous analyses (in this example called biomfile, and changed to biom_core) and fix metadata if needed (see Taxonomy_and_AlphaBetaDiv.md step 3.5 for instructions on this). 

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

**Step 2 - aggregate at taxonomic levels**

Make sure data is **compositional**, and aggregate at different taxonomic levels for different analyses. You can change “Phylum” to any other level of interest

```
core_rel.phy <- aggregate_taxa(core_rel, level="Phylum")
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

### Part 4 - ANCOM-BC calculations

Convert variables of interest to factors
```
ps.filt <- biomfile

sample_data(ps.filt)$Geographic_Location <- factor(
  sample_data(ps.filt)$Geographic_Location
)
sample_data(ps.filt)$Household_Location <- factor(
  sample_data(ps.filt)$Household_Location
)
```

Run ANCOM-BC analysis - change tax_level, fix_formula arguments, and group to fit your data

```
res <- ancombc2(
  data = ps.filt,
  tax_level = "Phylum",
  fix_formula = "Geographic_Location + Household_Location",
  group = "Geographic_Location",
  global = TRUE,
  pairwise = TRUE,
  dunnet = FALSE,
  trend = FALSE,
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  alpha = 0.05
)

# confirm data structure
names(res$res)

# generate tsv containing significant global results
global <- res$res_global
sig.global <- global %>%
  filter(q_val < 0.05)

write.table(
  sig.global,
  file = "ANCOMBC2_global_results_phylum_significant.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# generate tsv containing pairwise results
pairwise <- res$res_pair

write.table(
  pairwise,
  file = "ANCOMBC2_pairwise_results_phylum_all.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# generate summary tables from pairwise results for easier interpretation - reference level = B, you will need to change A, B, and C to match your variables. This code will be different with a larger number of variables.

summary_table <- pairwise %>%
  select(
    taxon,
    lfc_Geographic_LocationA,
    lfc_Geographic_LocationC,
    lfc_Geographic_LocationC_Geographic_LocationA
  ) %>%
  rename(
    Taxon = taxon,
    `A vs B` = lfc_Geographic_LocationA,
    `C vs B` = lfc_Geographic_LocationC,
    `C vs A` =
      lfc_Geographic_LocationC_Geographic_LocationA
  )

summary_table <- pairwise %>%
  left_join(
    global %>% select(taxon, q_val),
    by = "taxon"
  ) %>%
  select(
    taxon,
    q_val,
    lfc_Geographic_LocationA,
    lfc_Geographic_LocationC,
    lfc_Geographic_LocationC_Geographic_LocationA
  ) %>%
  rename(
    Taxon = taxon,
    `Global q-value` = q_val,
    `A vs B` = lfc_Geographic_LocationA,
    `C vs B` = lfc_Geographic_LocationC,
    `C vs A` =
      lfc_Geographic_LocationC_Geographic_LocationA
  )

sig_taxa <- global %>%
  filter(q_val < 0.05) %>%
  pull(taxon)

summary_table <- summary_table %>%
  filter(Taxon %in% sig_taxa)

# export tsv
write.table(
  summary_table,
  "ANCOMBC2_Phylum_summary_geoloc.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```
