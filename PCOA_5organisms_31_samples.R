# Load required libraries
library(vegan)      # For Bray-Curtis and PCoA
library(ggplot2)    # For plotting
library(dplyr)      # For data manipulation

# Step 1: Load abundance data and metadata
abundance_data <- read.table("Absolute_Abundance.txt", header = TRUE, row.names = 1, sep = "\t")
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t")

# Step 2: Align metadata with column names in abundance data
# Ensure organism names in metadata match column names in abundance data
metadata$Organism <- trimws(metadata$Organism)
colnames(abundance_data) <- trimws(colnames(abundance_data))

# Split organisms into Baseline and Treated groups and subset abundance data
baseline_data <- abundance_data[, colnames(abundance_data) %in% metadata$Organism[metadata$Group == "Baseline"]]
treated_data <- abundance_data[, colnames(abundance_data) %in% metadata$Organism[metadata$Group == "Treated"]]

# Load required libraries
library(vegan)      # For Bray-Curtis and PCoA
library(ggplot2)    # For plotting

# Step 1: Calculate Bray-Curtis dissimilarity matrices for Baseline and Treated
diss_baseline <- vegdist(baseline_data, method = "bray")
diss_treated <- vegdist(treated_data, method = "bray")

# Step 2: Perform PCoA on the dissimilarity matrices
pcoa_baseline <- cmdscale(diss_baseline, k = 2)
pcoa_treated <- cmdscale(diss_treated, k = 2)

# Step 3: Create data frames for plotting
pcoa_baseline_df <- data.frame(PC1 = pcoa_baseline[, 1], PC2 = pcoa_baseline[, 2], Group = "Baseline")
pcoa_treated_df <- data.frame(PC1 = pcoa_treated[, 1], PC2 = pcoa_treated[, 2], Group = "Treated")

# Step 4: Combine both PCoA results into one data frame
pcoa_df <- rbind(pcoa_baseline_df, pcoa_treated_df)

# Step 5: Plot the PCoA results
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  labs(title = "PCoA Plot (Bray-Curtis)", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = c("Baseline" = "blue", "Treated" = "red"))

