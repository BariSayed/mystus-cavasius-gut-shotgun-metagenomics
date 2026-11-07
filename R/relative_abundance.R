
# Taxonomic classification
# ====== Install required packages ======
install.packages(c("ggplot2", "dplyr", "tidyverse", "RColorBrewer", "patchwork"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("biomformat", "phyloseq"))

# ====== Install required packages (if needed) ======
# install.packages(c("ggplot2", "dplyr", "tidyverse", "RColorBrewer", "patchwork"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("biomformat", "phyloseq"))

# ====== Load Libraries ======
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(scales)
library(patchwork)

# ====== Import Data ======
physeq_raw <- import_biom("output.biom")
metadata <- read.csv("metadata.csv", row.names = 1, check.names = FALSE)
sample_data(physeq_raw) <- sample_data(metadata)

# ====== Clean Taxonomy ======
tax_table(physeq_raw) <- apply(tax_table(physeq_raw), 2, function(x) sub("^[a-z]__", "", x))
colnames(tax_table(physeq_raw)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ====== Filter Bacteria ======
physeq_filtered <- subset_taxa(physeq_raw, Kingdom == "Bacteria" & !grepl("Homo sapiens", Species, ignore.case = TRUE))
physeq_rel <- transform_sample_counts(physeq_filtered, function(x) x / sum(x))

# ====== Define swapped color palette ======
top15_colors <- c(
  "#2b8cbe",  # swapped: was #045a8d
  "#1b9e77",  # swapped: was #d7301f
  "#7bccc4",
  "#e34a33",
  "#045a8d",  # swapped: was #2b8cbe
  "#fec44f",
  "#fdbb84",
  "#fc8d59",
  "#78c679",
  "#41ab5d",
  "#d7301f",  # swapped: was #1b9e77
  "#fdae61",
  "#66c2a5",
  "#a6761d",
  "#666666"
)


# ====== PHYLUM ======
physeq_phylum <- tax_glom(physeq_rel, taxrank = "Phylum")
df_phylum <- psmelt(physeq_phylum)
df_phylum$Phylum <- as.character(df_phylum$Phylum)

top15_phyla <- df_phylum %>%
  group_by(Phylum) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:15) %>%
  pull(Phylum)

df_phylum <- df_phylum %>%
  mutate(Phylum = ifelse(Phylum %in% top15_phyla, Phylum, "Others"))

df_phylum_sum <- df_phylum %>%
  group_by(Sample, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance) * 100) %>%
  ungroup() %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

phylum_colors <- c(top15_colors, "gray70")
names(phylum_colors) <- c(top15_phyla, "Others")

phylum_plot <- ggplot(df_phylum_sum, aes(x = `Sample Name`, y = RelAbundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = phylum_colors) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(x = "Sample Replicates", y = "Relative Abundance (%)", title = "A", fill = "Phylum") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.text = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )


# ====== CLASS ======
# ====== CLASS ======
physeq_class <- tax_glom(physeq_rel, taxrank = "Class")
df_class <- psmelt(physeq_class)
df_class$Class <- as.character(df_class$Class)

top15_class <- df_class %>%
  group_by(Class) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:15) %>%
  pull(Class)

df_class <- df_class %>%
  mutate(Class = ifelse(Class %in% top15_class, Class, "Others"))

df_class_sum <- df_class %>%
  group_by(Sample, Class) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance) * 100) %>%
  ungroup() %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

class_colors <- c(top15_colors, "gray70")
names(class_colors) <- c(top15_class, "Others")

class_plot <- ggplot(df_class_sum, aes(x = `Sample Name`, y = RelAbundance, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(x = "Sample Replicates", y = "Relative Abundance (%)", title = "B", fill = "Class") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.text = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )

# ====== GENUS ======
# ====== GENUS ======
physeq_genus <- tax_glom(physeq_rel, taxrank = "Genus")
df_genus <- psmelt(physeq_genus)
df_genus$Genus <- as.character(df_genus$Genus)

top15_genus <- df_genus %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:15) %>%
  pull(Genus)

df_genus <- df_genus %>%
  mutate(Genus = ifelse(Genus %in% top15_genus, Genus, "Others"))

df_genus_sum <- df_genus %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance) * 100) %>%
  ungroup() %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

genus_colors <- c(top15_colors, "gray70")
names(genus_colors) <- c(top15_genus, "Others")

genus_plot <- ggplot(df_genus_sum, aes(x = `Sample Name`, y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = genus_colors) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(x = "Sample Replicates", y = "Relative Abundance (%)", title = "C", fill = "Genus") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13,, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 13, face = "bold")
  )


# ====== SPECIES ======
# ====== SPECIES ======
physeq_species <- tax_glom(physeq_rel, taxrank = "Species")
df_species <- psmelt(physeq_species)
df_species$GenusSpecies <- paste(df_species$Genus, df_species$Species)
df_species$GenusSpecies <- gsub("NA NA", "Unclassified", df_species$GenusSpecies)
df_species$GenusSpecies <- gsub(" NA", "", df_species$GenusSpecies)

top15_gs <- df_species %>%
  group_by(GenusSpecies) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>%
  slice(1:15) %>%
  pull(GenusSpecies)

df_species <- df_species %>%
  mutate(GenusSpecies = ifelse(GenusSpecies %in% top15_gs, GenusSpecies, "Others"))

df_species_sum <- df_species %>%
  group_by(Sample, GenusSpecies) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Abundance / sum(Abundance) * 100) %>%
  ungroup() %>%
  left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")

species_colors <- c(top15_colors, "gray70")
names(species_colors) <- c(top15_gs, "Others")

species_plot <- ggplot(df_species_sum, aes(x = `Sample Name`, y = RelAbundance, fill = GenusSpecies)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = species_colors) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(x = "Sample Replicates", y = "Relative Abundance (%)", title = "D", fill = "Genus Species") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 13, face = "bold")
  )

# ====== Combine Plots ======
combined_plot <- (phylum_plot | class_plot) / (genus_plot | species_plot)
combined_plot +
  plot_annotation(title = "Taxonomic Composition",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))

# Save high-resolution figure to your desired directory with increased size
ggsave(
  filename = "Taxonomic_Composition_Combined.png",
  plot = combined_plot,
  width = 18,      # increased width in inches
  height = 12,     # increased height in inches
  dpi = 800,       # high resolution
  units = "in"
)

