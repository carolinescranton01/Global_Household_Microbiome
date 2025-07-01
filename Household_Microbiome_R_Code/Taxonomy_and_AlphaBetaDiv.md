## RMarkdown walkthrough on how to analyze .biom files for taxonomic composition, alpha diversity, and beta diversity

This code was used to analyze all samples for both bacterial, viral, and eukaryotic pathogen taxa, with a few minor edits (such as changing the taxonomic ranks for viruses)

### Part 1 - setup and data cleanup

**Step 1. Load required packages**

```
library(speedyseq)
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

# Alternativly, fulter to only include kingdom Bacteria
biomfile <- subset_taxa(biomfile, Kingdom == "Bacteria")
```

Now we should just be left with bacteria reads, and can start analyzing these reads

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
```

**NOTE**: Alpha diversity was exported into excel and added into a new dataframe with sample metadata (sample ID, geographic location, household location, etc) which was re-uploaded to R and used for ggplot box+whisker plots which used ANOVA for statistical analysis. 

Below is an example of the code for a plot with statistics. X can be changed to different columns in the metadata, and Y can be changed to different alpha diversity metrics, calculated with the command above:

```
household_shannon_div <- ggboxplot(diversity_dataframe, x = "Household_Location", y = "diversity_shannon") +
    rotate_x_text() +
    ylim(0, 10) +
    theme(legend.position="none") +
    labs(title="Shannon Diversity in Households in Different Cities") +
    stat_compare_means(aes(group=Household_Location), label = 'p.format', method='anova', label.y=5, label.x=1.5) +
    facet_wrap(~Geographic_Location)
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
biom_phylum.rel <- psmelt(biom_phylum.rel)

# Convert phyloseq object to a data frame
dfphy <- psmelt(biom_phylum.rel)
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
