## Antibiotic resistance gene and virulence factor abundance and alpha/beta diversity analysis

## Part 1 - Setup

**Step 1 - load required packages**

```
library(ggplot2) # Version 4.0.2
library(ggpubr) # Version 0.6.3
library(tidyverse) # Version 2.0.0
library(broom) # Version 1.0.12
library(AICcmodavg) # Version 2.3.4
library(readxl) # Version 1.4.5
library(rstatix) # Version 1.7.3
library(microbiome) # Version 1.30.0
library(dplyr) # Version 1.2.0
library(writexl) # Version 1.5.4
library(multcompView) # Version 0.1.11
```

**Step 2 - re-format abricate outputs (.tab files, one per sample per database) into excel sheets - outside of R**

One sheet was made for the virulence factors (detected using the Virulence Factor Database VFDB), and another for the antibiotic resistance genes detected with the ResFinder database. Columns included the Sample_ID and other metadata (location, house number, household location, scaled number of reads) and then gene name and gene function as two columns, followed by the number of occurrences of that gene in the specific sample, and the scaled occurrences (occurrences / scaled reads) in the sample. Additional columns containing metadata for the samples were added by hand - these columns included household number, household location, geographic location, number of reads (total), scaled reads (total reads/1000), and scaled occurrences of each gene (occurrences of gene / scaled reads). All of this reformatting was done in excel, using the python scripts writeexcel.py and sortexcel.py to combine the data (found in this github repository), metadata was added by hand, and then the consolidated_data.xslx files with metadata were imported into R for analysis using the read_excel function. Additionally, short descriptions were added as a column in this table for graphing.

## Part 2 - Abundance Analysis

**Step 1 - subset out the top 20 most-abundant genes**

The top 20 VFs and ARGs were used in this analysis, however you could look at the top 5, top 10, top 100, etc - just change '20' to the number of choice
The example dataframe is called VF, with sheet_name (sample ID) gene_graph (gene name and function info for plot) and scaled_occurrences (calculated in excel by dividing the number of occurances by the scaled number of reads in each sample) columns. Any dataframe with columns to denote the sample ID, gene name, gene function, and scaled occurrences can be used.

```
# rename columns for plotting
VF_plot <- VF %>%
  rename(
    Sample = sheet_name,
    Gene = gene_graph,
    Abundance = scaled_occurances
  )

# top 20
top20_VF <- VF_plot %>%
  group_by(Gene) %>%
  summarise(total = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  slice_head(n = 20) %>%
  pull(Gene)

# create other category
top20_VF <- top20_VF %>%
  mutate(Gene_grouped = ifelse(Gene %in% top20_VF, Gene, "Other"))

# rename column to sampleID for metadata joining
colnames(top20_VF)[1] <- "Sample_ID"

VF_df_annotated <- top20_VF %>%
  left_join(metadata, by = "Sample_ID")

```

**Step 2 - Generating relative abundance data for top 20 VFs/ARGs**

Factors in the group_by() argument are different metadata variables - in this example, this grouping was structured so that samples were all looked at individually - if more than one type of gene was found in a sample, these two or more genes from the same sample needed to be looked at together to determine their relative abundance. The grouping structure below is Geographic_Location > Household_Number > Household_Location, so samples from location 1 in house 1 of country 1 would be grouped together, samples from location 2 in house 1 in country 1 would be together, and so on.

```
# prepare data for plot by grouping variables for relative abundance calculations
VF_df_annotated %>%
  group_by(Geographic_Location, Household_Location, Gene) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE),
            .groups = "drop")

# calculate relative abundance
VF_df_rel <- VF_df_annotated %>%
  group_by(Geographic_Location, Household_Location) %>%
  mutate(RelAbundance = Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

# keep only complete cases
VF_df_rel <- VF_df_rel[complete.cases(VF_df_rel), ]

```

**Step 4 - make a relative abundance bar graph, showing the relative abundance (out of 1.0) of the top 20 genes in each sample**

To make a raw abundance graph, replace y = RelAbundance with scaled_occurrences

```
# plot with ggplot - change any relevant variables
plot <- ggplot(VF_df_rel, aes(
  x = Household_Location,
  y = RelAbundance,
  fill = Gene_grouped)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.2, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    strip.text = element_text(size = 12)
  ) +
  labs(
    title = "Relative Abundance of Virulence Factor Genes Found in Different Locations",
    x = "Household Location",
    y = "Relative Abundance",
    fill = "Gene"
  ) +
  facet_wrap(~Geographic_Location) +
  scale_x_discrete(limits = unique(res_df_rel$Household_Location)) +
  scale_fill_manual(values = kelly_colors) +
  guides(fill = guide_legend(ncol = 4)
)

print(plot)
```

## Part 3 - Gene Alpha and Beta Diversity

**Step 1 - data restructuring**

This requires certain data to be restructured again, but in R this time. Below is an example on how to restructure the data. The example data is called biom_res and is a data frame with 3 columns – FullName (which is the sample_ID), Gene_Name, and Scaled_Occurrences. If one sample had ten different ARGs/VFs detected, it will have ten rows in this data frame with the same FullName value and different Gene_Name and Scaled_Occurrences values for each row. Data was structured in Microsoft Excel and imported using read_xlsx(). This also requires a metadata file with info on the samples to use when graphing diversity (ie. sample ID, location) – in this example called metadata_res.

Restructure data so that there are now many columns (each is a Gene_Name), and all samples (FullName) have one row. Should be a large matrix, where if a specific gene is detected in a sample, it’s scaled occurrence will be in that column, but if a gene was not detected in that sample, it will have a 0. 

Data in this example is from Resfinder, in a dataframe called biom_res, and has columns sheet_name (sample ID), Gene_Name (description of gene for graph), and Scaled_Occurrences with metadata in a dataframe called metadata_res.

<img width="374" alt="Screenshot 2025-05-12 at 3 38 16 PM" src="https://github.com/user-attachments/assets/6e46db5a-289f-47a6-be69-729eb2a9838f" />

Run the pivot_wider command:

```
wide_data_res <- biom_res %>%
        pivot_wider(names_from = Gene_Name, values_from = Scaled_Occurrences, values_fill=0)
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
    left_join(metadata_res, by = "SampleID")
merged_results_res <- na.omit(merged_results_res)
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

# Tukey's post-hoc test, to see which variables significantly differ:
geo_aov <- aov(shannon_diversity ~ Geographic_Location, data = merged_results_res)
tuk <- TukeyHSD(arg_geo_aov)
letters <- multcompLetters4(arg_geo_aov, tuk)
print(letters)
```

**Step 4 - Beta diversity analysis**

Beta diversity on genes uses the same re-formatted data as above. In this example data is called wide_data_res (same as the alpha diversity's wide data). Metadata is structured the same, again only for composites in this example, as only composite data was used for this particular analysis. 

Generating bray-curtis distances, formatting as a data frame, generate PC values for axis, and attach to metadata: 

```
bray_curtis_dist <- vegdist(wide_data_res[,-1], method = "bray")
dist_matrix <- as.matrix(bray_curtis_dist)

pcoa_results_resbeta <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)

pcoa_df_resbeta <- data.frame(PC1 = pcoa_results_resbeta$points[,1], PC2 = pcoa_results_resbeta$points[,2], Sample_ID = wide_data_res$sheet_name)
colnames(pcoa_df_resbeta)[3] <- "Sample_ID"
pcoa_df_resbeta <- pcoa_df_resbeta %>%
  left_join(metadata, by = "Sample_ID")

pc1_perc <- round(
  pcoa_results_resbeta$eig[1] / sum(pcoa_results_resbeta$eig) * 100,
  1
)

pc2_perc <- round(
  pcoa_results_resbeta$eig[2] / sum(pcoa_results_resbeta$eig) * 100,
  1
)
```

Example code for PCoA plot:

```
res_beta <- ggplot(
  pcoa_df_resbeta,
  aes(x = PC1, y = PC2, color = Geographic_Location)) +
  geom_point(size = 3) +
  labs(
    title = "PCoA of ARGs in Global Cities",
    x = paste0("PC1 (", pc1_perc, "%)"),
    y = paste0("PC2 (", pc2_perc, "%)"),
    color = "Location" ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
  ) +
  scale_color_manual(values = pal)

print(res_beta)
```

Running a PERMANOVA test on wide data and test for beta dispersion via ANOVA

```
adonis1 <- adonis2(
  dist_matrix ~ Household_Location,
  data = metadata_res,
  by = "margin"
)

dist_matrix_dist <- as.dist(dist_matrix)

disp <- betadisper(
  dist_matrix_dist,
  metadata_res$Household_Location
)

print(adonis1)
anova(disp)
```













