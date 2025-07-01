## Antibiotic resistance gene and virulence factor abundance and alpha/beta diversity analysis

## Part 1 - Setup

**Step 1 - load required packages**

```
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(readxl)
library(rstatix)
library(microbiome)
library(dplyr)
library(writexl)
```

**Step 2 - re-format abricate outputs (.tab files, one per sample per database) into excel sheets - outside of R**

One sheet was made for the virulence factors (detected using the Virulence Factor Database VFDB), and another for the antibiotic resistance genes detected with the ResFinder database. Columns included the Sample_ID and other metadata (location, house number, household location, scaled number of reads) and then gene name and gene function as two columns, followed by the number of occurrences of that gene in the specific sample, and the scaled occurrences (occurrences / scaled reads) in the sample. Additional columns containing metadata for the samples were added by hand - these columns included household number, household location, geographic location, number of reads (total), scaled reads (total reads/1000), and scaled occurrences of each gene (occurrences of gene / scaled reads). All of this reformatting was done in excel, using the python scripts writeexcel.py and sortexcel.py to combine the data (found in this github repository), metadata was added by hand, and then the consolidated_data.xslx files with metadata were imported into R for analysis using the read_excel function. 

## Part 2 - Abundance Analysis

**Step 1 - subset out the top 20 most-abundant genes**

The top 20 VFs and ARGs were used in this analysis, however you could look at the top 5, top 10, top 100, etc - just change '20' to the number of choice
The example dataframe is called VF, with gene_function and scaled_occurrences columns. Any dataframe with columns to denote the sample ID, gene name, gene function, and scaled occurrences can be used.

```
virulence_20 <- VF %>%
    group_by(gene_function) %>%
    summarise(scaled_occurrences = sum(scaled_occurrences, na.rm = TRUE)) %>%
    arrange(desc(scaled_occurrences))
top_20_factors_VF <- virulence_20 %>%
    top_n(20, scaled_occurrences)
filtered_VFtop20 <- VF %>%
    filter(gene_function %in% top_20_factors_VF$gene_function)
print(top_20_factors_VF)
print(filtered_VFtop20)
```

**Step 2 - Generating relative abundance data for top 20 VFs/ARGs**

Factors in the group_by() argument are different metadata variables - in this example, this grouping was structured so that samples were all looked at individually - if more than one type of gene was found in a sample, these two or more genes from the same sample needed to be looked at together to determine their relative abundance. The grouping structure below is Geographic_Location > Household_Number > Household_Location, so samples from location 1 in house 1 of country 1 would be grouped together, samples from location 2 in house 1 in country 1 would be together, and so on.

```
top_relative_abundance_VFDB_comp <- filtered_VFDBtop20 %>%
    group_by(Geographic_Location, Household_Number, Household_Location) %>%
    mutate(total_occurrences = sum(occurrences)) %>%
    mutate(Relative_Abundance = occurrences / total_occurrences) %>%
    ungroup() %>%
    select(Geographic_Location, Household_Number, Household_Location, Gene_Name,
    Gene_function, scaled_occurrences, Relative_Abundance)
```

**Step 3 - Summarize and filter data to only keep data which exists in the processed set (the top 20 genes and their relative abundances)**

```
summarized_data <- top_relative_abundance_VFDB_comp %>%
    group_by(Geographic_Location, Household_Location, Gene_function) %>%
    summarise(Total_Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE)) %>%
    ungroup()
complete_data <- summarized_data %>%
    complete(Geographic_Location, Household_Location)
    complete_data[is.na(complete_data$Total_Relative_Abundance), "Total_Relative_Abundance"] <- 0
filtered_data <- complete_data %>% filter(Total_Relative_Abundance > 0)
```

**Step 4 - make a relative abundance bar graph, showing the relative abundance (out of 1.0) of the top 20 genes in each sample**

To make an actual abundance graph, replace y = Total_Relative_Abundance with scaled_occurrences


```
VFDB_RA_barplot_comp <- ggplot(filtered_data, aes(x = Household_Location, y = Total_Relative_Abundance, fill = Gene_function)) + 
    geom_bar(stat = "identity", position = "stack", na.rm = TRUE, alpha = 1, color = "black", linewidth = 0.3, show.legend = TRUE) +
  	font("ylab", size = 15, face = "bold") +
 	  font("xlab", size = 15, face = "bold") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(title = "Relative Abundance of VF Genes Found in Different Locations", size = 30) +
  	ylim(0, 1) + 
  	theme_bw() + 
  	facet_wrap(~Geographic_Location) + 
    scale_x_discrete(limits = unique(filtered_data$Household_Location)) +
  	theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 15, face = "bold"),   axis.title.y = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), legend.position = "bottom", legend.title = element_text(size = 12),
        legend.text = element_text(size = 10), plot.margin = margin(1, 1, 2, 1, "cm"), strip.text = element_text(size = 12, face = "bold"))

print(VFDB_RA_barplot_comp)
```

## Part 3 - Gene Alpha and Beta Diversity

**Step 1 - data restructuring**

This requires certain data to be restructured again, but in R this time. Below is an example on how to restructure the data. The example data is called biom_res and is a data frame with 3 columns – FullName (which is the sample_ID), Gene_Name, and Scaled_Occurrences. If one sample had ten different ARGs/VFs detected, it will have ten rows in this data frame with the same FullName value and different Gene_Name and Scaled_Occurrences values for each row. Data was structured in Microsoft Excel and imported using read_xlsx(). This also requires a metadata file with info on the samples to use when graphing diversity (ie. sample ID, location) – in this example called metadata_res.

Restructure data so that there are now many columns (each is a Gene_Name), and all samples (FullName) have one row. Should be a large matrix, where if a specific gene is detected in a sample, it’s scaled occurrence will be in that column, but if a gene was not detected in that sample, it will have a 0. 

<img width="374" alt="Screenshot 2025-05-12 at 3 38 16 PM" src="https://github.com/user-attachments/assets/6e46db5a-289f-47a6-be69-729eb2a9838f" />

Run the pivot_wider command:

```
wide_data_res <- biom_res %>%
        pivot_wider(names_from = Gene_Name, values_from = Scaled_Occurances, values_fill=0)
```
<img width="376" alt="Screenshot 2025-05-12 at 3 38 26 PM" src="https://github.com/user-attachments/assets/c8f47246-49a2-44e9-bf1b-2e86ce6acca2" />


**Step 2 - alpha diversity**

Run alpha diversity calculations on the restructured data, convert it to data frame, and merge it with the metadata (FullName is the shared sample ID between the restructed data and the metadata dataframe)

```
diversity_results_res <- wide_data_res %>%
    select(-FullName) %>%
    apply(1, function(x) diversity(x, index = "shannon"))

results_res <- data.frame(FullName = wide_data_res$FullName,
    shannon_diversity=diversity_results)
merged_results_res <- results_res %>%
    left_join(metadata_res, by = "FullName")
```

**Step 3 - box-and-whisker plots of alpha diversity metrics, with ANOVA**

Example ggplot code to make box and whisker plots with ANOVA statistical analysis using the data generated above. Change variables as needed to fit your data/metadata/alpha diversity metric

```
res_geolocs_alpha <- ggboxplot(merged_results_res, x = "Geographic_Location", y = "shannon_diversity") +
    rotate_x_text() +
    ylim(0, 2) +
    theme(legend.position="none") +
    labs(title="Shannon Diversity in ARGs in Different Cities") +
    stat_compare_means(aes(group=Geographic_Location), label = 'p.format', method='anova', label.y=1.99, label.x=0.6)
print(res_geolocs_alpha)
```

**Step 4 - Beta diversity analysis**

Beta diversity on genes uses the same re-formatted data as above. In this example data is called wide_data_res_beta – similar to wide_data_res, but in this example used **composite data** (combined samples from each household location per geographic location, ie all coffee tables in Tucson, AZ) rather than the entire non-composite set. Metadata is structured the same, again only for composites in this example, as only composite data was used for this particular analysis. 

Generating bray-curtis distances, formatting as a data frame, and linking to metadata: 

```
bray_curtis_dist_resbeta <- vegdist(wide_data_res_beta[,-1], method = "bray")
pcoa_results_resbeta <- cmdscale(bray_curtis_dist_resbeta, eig = TRUE, k = 2)
pcoa_df_resbeta <- data.frame(PC1 = pcoa_results_resbeta$points[,1], PC2 = >pcoa_results_resbeta$points[,2], FullName = wide_data_res_beta$FullName)
pcoa_df_resbeta <- pcoa_df_resbeta %>%
        left_join(metadata_res_beta, by = "FullName")
```

Example code for PCoA plot:

```
res_beta <- ggplot(pcoa_df_resbeta, aes(x = PC1, y = PC2, color = Household_Location, shape=Geographic_Location)) +
        geom_point(size = 3) +
        labs(title = "PCoA of Bray-Curtis Dissimilarity", x = "PC1", y = "PC2") +
        theme_minimal() +
        stat_ellipse(aes(color = Geographic_Location), level = 0.95, size = 0.5)
print(res_beta)
```

Running a PERMANOVA test with wide_data and metadata (matched sampleIDs with any metadata of interest):

```
# Remove any NAs and replace with 0
wide_data_res <- wide_data_res %>%
  filter(rowSums(is.na(.)) == 0)

# Calculate bray curtis distances (not counting first column, which is sampleIDs)
bray_curtis_dist <- vegdist(wide_data_res[,-1], method = "bray")

# Convert distances to matrix
dist_matrix <- as.matrix(bray_curtis_dist)

# Re-add the sampleIDs to the distance matrix
rownames(dist_matrix) <- wide_data_res$SampleID

# Verify the two items have the same sampleIDs
cat("Samples in distance matrix:", head(rownames(dist_matrix)), "\n")
cat("Samples in metadata:", head(metadata), "\n")

# Verify that the two items have the same number of rows
cat("Rows in distance matrix:", nrow(dist_matrix), "\n")
cat("Rows in metadata:", nrow(metadata), "\n")

# Convert matrix to distance for PERMANOVA
dist <- as.dist(dist_matrix)

# Conduct PERMANOVA test (change Household_Location to metadata variable of interest)
PERMANOVA <- adonis2(dist ~ Household_Location, 
                       data = metadata)
print(PERMANOVA)
```













