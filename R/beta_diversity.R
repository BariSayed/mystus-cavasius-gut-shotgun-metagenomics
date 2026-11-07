# Install packages
install.packages("vegan")

# Load libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)

# Set working directory
# setwd(".")  # set working directory to the folder containing output.biom and metadata.csv

# Import biom file and metadata
merged_metagenomes <- import_biom("output.biom")

# Clean taxonomy labels (remove prefixes)
tax_table(merged_metagenomes) <- apply(tax_table(merged_metagenomes), 2, function(x) sub("^[a-z]__", "", x))
colnames(tax_table(merged_metagenomes)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Subset only bacteria
merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria")

# Load metadata and clean
metadata <- read.csv("metadata.csv", row.names = 1)
rownames(metadata) <- trimws(rownames(metadata))

# Attach metadata to phyloseq object
sample_data(merged_metagenomes) <- sample_data(metadata)

# Calculate Bray-Curtis distance matrix (common for beta diversity)
bray_dist <- phyloseq::distance(merged_metagenomes, method = "bray")

# PERMANOVA test to check group effect on beta diversity
adonis_res <- vegan::adonis2(bray_dist ~ Group, data = metadata)
print(adonis_res)

# ---- PCoA Plot ----
pcoa_res <- ape::pcoa(bray_dist)

# Extract PCoA scores for plotting
pcoa_df <- data.frame(pcoa_res$vectors[,1:2], SampleID = rownames(pcoa_res$vectors))
pcoa_df <- left_join(pcoa_df, metadata %>% tibble::rownames_to_column("SampleID"), by = "SampleID")

# % variance explained for axis labels
var_exp <- round(pcoa_res$values$Relative_eig * 100, 1)

pcoa_plot <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size = 4) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.3, color = NA) +
  labs(
    title = "A. PCoA based on Bray-Curtis",
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ---- NMDS Plot ----
set.seed(123)  # for reproducibility
nmds_res <- metaMDS(bray_dist, k = 2, trymax = 100)

nmds_df <- data.frame(NMDS1 = nmds_res$points[,1], NMDS2 = nmds_res$points[,2], SampleID = rownames(nmds_res$points))
nmds_df <- left_join(nmds_df, metadata %>% tibble::rownames_to_column("SampleID"), by = "SampleID")

stress <- round(nmds_res$stress, 4)

nmds_plot <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Group)) +
  geom_point(size = 4) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.3, color = NA) +
  labs(
    title = "B. NMDS based on Bray-Curtis",
    x = "NMDS1",
    y = "NMDS2",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# Combine PCoA and NMDS plots side-by-side
combined_beta_plot <- pcoa_plot | nmds_plot

print(combined_beta_plot)

# Save combined figure as high-res PNG
ggsave("Beta_Diversity_Combined_PCoA_NMDS.png", combined_beta_plot, width = 12, height = 6, dpi = 600)


var_exp <- round(pcoa_res$values$Relative_eig * 100, 1)
print(var_exp[1:2])  # PCoA1 and PCoA2

stress <- round(nmds_res$stress, 4)
print(stress)

adonis_res <- vegan::adonis2(bray_dist ~ Group, data = metadata)
print(adonis_res)

adonis_R2 <- adonis_res$R2[1]
adonis_p <- adonis_res$`Pr(>F)`[1]
print(paste("R² =", round(adonis_R2, 3), "; p =", adonis_p))


beta_disp <- betadisper(bray_dist, metadata$Group)
permutest_res <- permutest(beta_disp)
print(permutest_res)
