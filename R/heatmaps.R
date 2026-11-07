
# ========== Load Libraries ==========
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)  # For Brewer palettes
library(tibble)
library(patchwork)

# ========== Set Working Directory ==========
# setwd(".")  # set working directory to the folder containing output.biom and metadata.csv

# ========== Load Metadata ==========
metadata <- read.csv("metadata.csv", row.names = 1, check.names = FALSE)

# ========== Load BIOM and Attach Metadata ==========
physeq <- import_biom("output.biom")
sample_data(physeq) <- sample_data(metadata)

# ========== Clean Taxonomy ==========
tax_table(physeq) <- apply(tax_table(physeq), 2, function(x) sub("^[a-z]__", "", x))
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ========== Subset Bacteria Only ==========
physeq_bact <- subset_taxa(physeq, Kingdom == "Bacteria")

# ========== Transform to Relative Abundance ==========
physeq_rel <- transform_sample_counts(physeq_bact, function(x) x / sum(x))

# ========== A. Top 50 Genus ==========
phy_genus <- tax_glom(physeq_rel, taxrank = "Genus")
df_genus <- psmelt(phy_genus)

top50_genus <- df_genus %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:50) %>%
  pull(Genus)

df_genus <- df_genus %>%
  filter(Genus %in% top50_genus) %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

df_genus$`Sample Name` <- factor(df_genus$`Sample Name`, levels = c("CR1", "CR2", "TR1", "TR2", "TR3"))

genus_heat <- ggplot(df_genus, aes(x = `Sample Name`, y = Genus, fill = Abundance)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_distiller(palette = "Paired", name = "Relative Abundance", direction = 1) +
  labs(title = "A. Top 50 Genera", x = "Conditions", y = "Genus") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0, margin = margin(b = 10)),
    axis.title.x = element_text(face = "bold", size = 18, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 18, margin = margin(r = 10)),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 15),
    panel.grid.major.y = element_line(color = "gray85", size = 0.3),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 13)
  )

# ========== B. Top 50 Species ==========
phy_species <- tax_glom(physeq_rel, taxrank = "Species")
df_species <- psmelt(phy_species)

# Create "Scientific" column as Genus + Species
df_species <- df_species %>%
  mutate(Scientific = paste(Genus, Species, sep = " ")) %>%
  group_by(Scientific) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:50) %>%
  pull(Scientific) -> top50_species

df_species <- psmelt(phy_species) %>%
  mutate(Scientific = paste(Genus, Species, sep = " ")) %>%
  filter(Scientific %in% top50_species) %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

df_species$`Sample Name` <- factor(df_species$`Sample Name`, levels = c("CR1", "CR2", "TR1", "TR2", "TR3"))

species_heat <- ggplot(df_species, aes(x = `Sample Name`, y = Scientific, fill = Abundance)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_distiller(palette = "Paired", name = "Relative Abundance", direction = 1) +
  labs(title = "B. Top 50 Species", x = "Conditions", y = "Scientific Name") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0, margin = margin(b = 10)),
    axis.title.x = element_text(face = "bold", size = 18, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 18, margin = margin(r = 10)),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 15),
    panel.grid.major.y = element_line(color = "gray85", size = 0.3),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 13)
  )

# ========== Combine Heatmaps Side-by-Side ==========
combined_plot <- genus_heat | species_heat + plot_layout(guides = "collect")

# ========== Print and Save ==========
print(combined_plot)

# Save
ggsave("Top50_Genus_Species_Heatmap_Paired.png", 
       plot = combined_plot, width = 18, height = 12, dpi = 600)

cat("✅ Heatmaps for top 50 genera and species saved as 'Top50_Genus_Species_Heatmap_Paired.png' in the specified folder.\n")









---------------------
  
# ========== Load Libraries ==========
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)  # For Brewer palettes
library(viridis)       # If you still want to compare
library(tibble)
library(patchwork)

# ========== Set Working Directory ==========
# setwd(".")  # set working directory to the folder containing output.biom and metadata.csv

# ========== Load Metadata ==========
metadata <- read.csv("metadata.csv", row.names = 1, check.names = FALSE)

# ========== Load BIOM and Attach Metadata ==========
physeq <- import_biom("output.biom")
sample_data(physeq) <- sample_data(metadata)

# ========== Clean Taxonomy ==========
tax_table(physeq) <- apply(tax_table(physeq), 2, function(x) sub("^[a-z]__", "", x))
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ========== Subset Bacteria Only ==========
physeq_bact <- subset_taxa(physeq, Kingdom == "Bacteria")

# ========== Transform to Relative Abundance ==========
physeq_rel <- transform_sample_counts(physeq_bact, function(x) x / sum(x))

# ========== A. Top 30 Phylum ==========
phy_phylum <- tax_glom(physeq_rel, taxrank = "Phylum")
df_phylum <- psmelt(phy_phylum)

top30_phylum <- df_phylum %>%
  group_by(Phylum) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:30) %>%
  pull(Phylum)

df_phylum <- df_phylum %>%
  filter(Phylum %in% top30_phylum) %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

df_phylum$`Sample Name` <- factor(df_phylum$`Sample Name`,
                                  levels = c("CR1", "CR2", "TR1", "TR2", "TR3"))

phylum_heat <- ggplot(df_phylum, aes(x = `Sample Name`, y = Phylum, fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "Paired", name = "Rel. Abundance", direction = 1) +
  labs(title = "a. Top 30 Phyla", x = "Conditions", y = "Phylum") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 16),              # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 16),              # Bold y-axis title
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),  # Bold + bigger CR1, CR2...
    axis.text.y = element_text(face = "bold", size = 10),               # Bold + bigger Phylum names
    plot.title = element_text(face = "bold", hjust = 0)
  )

# ========== B. Top 30 Genus ==========
phy_genus <- tax_glom(physeq_rel, taxrank = "Genus")
df_genus <- psmelt(phy_genus)

top30_genus <- df_genus %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:30) %>%
  pull(Genus)

df_genus <- df_genus %>%
  filter(Genus %in% top30_genus) %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

df_genus$`Sample Name` <- factor(df_genus$`Sample Name`,
                                 levels = c("CR1", "CR2", "TR1", "TR2", "TR3"))

genus_heat <- ggplot(df_genus, aes(x = `Sample Name`, y = Genus, fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "Paired", name = "Rel. Abundance", direction = 1) +
  labs(title = "b. Top 30 Genera", x = "Conditions", y = "Genus") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 16),              # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 16),              # Bold y-axis title
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),  # Bold + bigger CR1, CR2...
    axis.text.y = element_text(face = "bold", size = 10),               # Bold + bigger Genus names
    plot.title = element_text(face = "bold", hjust = 0)
  )

# ========== Combine Heatmaps Side-by-Side ==========
combined_plot <- phylum_heat | genus_heat + plot_layout(guides = "collect")

# ========== Print and Save ==========
print(combined_plot)
ggsave("combined_heatmaps_paired.png", combined_plot, width = 12, height = 8, dpi = 600)













--------------------------------------------------------
# Load libraries
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(tibble)
library(patchwork)

# Set working directory
# setwd(".")  # set working directory to the folder containing output.biom and metadata.csv

# Load metadata
metadata <- read.csv("metadata.csv", row.names = 1, check.names = FALSE)

# Load biom file
physeq <- import_biom("output.biom")
sample_data(physeq) <- sample_data(metadata)

# Clean taxonomy
tax_table(physeq) <- apply(tax_table(physeq), 2, function(x) sub("^[a-z]__", "", x))
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Subset to Bacteria only
physeq_bact <- subset_taxa(physeq, Kingdom == "Bacteria")

# Relative abundance
physeq_rel <- transform_sample_counts(physeq_bact, function(x) x / sum(x))

# ========== A. Top 30 Phylum ==========
phy_phylum <- tax_glom(physeq_rel, taxrank = "Phylum")
df_phylum <- psmelt(phy_phylum)

top30_phylum <- df_phylum %>%
  group_by(Phylum) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:30) %>%
  pull(Phylum)

df_phylum <- df_phylum %>%
  filter(Phylum %in% top30_phylum) %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

df_phylum$`Sample Name` <- factor(df_phylum$`Sample Name`,
                                  levels = c("CR1", "CR2", "TR1", "TR2", "TR3"))

phylum_heat <- ggplot(df_phylum, aes(x = `Sample Name`, y = Phylum, fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_viridis(name = "Rel. Abundance", option = "C") +
  labs(title = "A. Top 30 Phyla", x = NULL, y = "Phylum") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 6),
    plot.title = element_text(face = "bold", hjust = 0)
  )

# ========== B. Top 30 Genus ==========
phy_genus <- tax_glom(physeq_rel, taxrank = "Genus")
df_genus <- psmelt(phy_genus)

top30_genus <- df_genus %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:30) %>%
  pull(Genus)

df_genus <- df_genus %>%
  filter(Genus %in% top30_genus) %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

df_genus$`Sample Name` <- factor(df_genus$`Sample Name`,
                                 levels = c("CR1", "CR2", "TR1", "TR2", "TR3"))

genus_heat <- ggplot(df_genus, aes(x = `Sample Name`, y = Genus, fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_viridis(name = "Rel. Abundance", option = "C") +
  labs(title = "B. Top 30 Genera", x = NULL, y = "Genus") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 6),
    plot.title = element_text(face = "bold", hjust = 0)
  )

# ========== Combine Side-by-Side ==========
combined_plot <- phylum_heat | genus_heat + plot_layout(guides = "collect")

# Print and save high-res figure
print(combined_plot)
ggsave("Top30_Phylum_Genus_Heatmap.png", plot = combined_plot, width = 14, height = 10, dpi = 600)

------------------
# Load libraries
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(patchwork)

# Load metadata
metadata <- read.csv("metadata.csv", row.names = 1, check.names = FALSE)

# Load biom file
physeq <- import_biom("output.biom")
sample_data(physeq) <- sample_data(metadata)

# Clean taxonomy
tax_table(physeq) <- apply(tax_table(physeq), 2, function(x) sub("^[a-z]__", "", x))
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Subset to Bacteria only
physeq_bact <- subset_taxa(physeq, Kingdom == "Bacteria")

# Relative abundance
physeq_rel <- transform_sample_counts(physeq_bact, function(x) x / sum(x))

# Define smart, interpretable color palette (blue to red)
smart_colors <- c("#e5f5f9", "#9ecae1", "#fdcc8a", "#fc8d59", "#d7301f")

# ========== A. Top 30 Phylum ==========
phy_phylum <- tax_glom(physeq_rel, taxrank = "Phylum")
df_phylum <- psmelt(phy_phylum)

top30_phylum <- df_phylum %>%
  group_by(Phylum) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:30) %>%
  pull(Phylum)

df_phylum <- df_phylum %>%
  filter(Phylum %in% top30_phylum) %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

df_phylum$`Sample Name` <- factor(df_phylum$`Sample Name`,
                                  levels = c("CR1", "CR2", "TR1", "TR2", "TR3"))

phylum_heat <- ggplot(df_phylum, aes(x = `Sample Name`, y = Phylum, fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = smart_colors, name = "Rel. Abundance") +
  labs(title = "A. Top 30 Phyla", x = NULL, y = "Phylum") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold", size = 12, hjust = 0)
  )

# ========== B. Top 30 Genus ==========
phy_genus <- tax_glom(physeq_rel, taxrank = "Genus")
df_genus <- psmelt(phy_genus)

top30_genus <- df_genus %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:30) %>%
  pull(Genus)

df_genus <- df_genus %>%
  filter(Genus %in% top30_genus) %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

df_genus$`Sample Name` <- factor(df_genus$`Sample Name`,
                                 levels = c("CR1", "CR2", "TR1", "TR2", "TR3"))

genus_heat <- ggplot(df_genus, aes(x = `Sample Name`, y = Genus, fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = smart_colors, name = "Rel. Abundance") +
  labs(title = "B. Top 30 Genera", x = NULL, y = "Genus") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold", size = 12, hjust = 0)
  )

# ========== Combine Side-by-Side ==========
combined_plot <- phylum_heat | genus_heat + plot_layout(guides = "collect")

# Print and save high-res figure
print(combined_plot)
ggsave("Top30_Phylum_Genus_Heatmap_SmartColor.png", plot = combined_plot, width = 14, height = 10, dpi = 600)

-----------------------------------------------------------------------
install.packages("viridis")  # Install
library("viridis") 
install.packages("RColorBrewer")
library(RColorBrewer)
display.brewer.all()


