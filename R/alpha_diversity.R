# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(FSA)          # For Dunn test
library(rcompanion)   # For compact letter display
library(ggpubr)

# Set working directory
# setwd(".")  # set working directory to the folder containing output.biom and metadata.csv

# Import BIOM file
physeq <- import_biom("output.biom")

# Clean taxonomic labels
tax_table(physeq) <- apply(tax_table(physeq), 2, function(x) sub("^[a-z]__", "", x))
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Filter for Bacteria only
physeq <- subset_taxa(physeq, Kingdom == "Bacteria")

# Estimate alpha diversity
alpha_df <- estimate_richness(physeq, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
alpha_df$SampleID <- rownames(alpha_df)
alpha_df$InvSimpson <- 1 / alpha_df$Simpson

# Load metadata with SampleID as a column
metadata <- read.csv("metadata.csv", check.names = FALSE)

# Merge alpha diversity with metadata by SampleID
alpha_df <- merge(alpha_df, metadata, by = "SampleID")

# Statistical significance test function
get_signif_letters <- function(df, metric, group_col = "Group") {
  kw <- kruskal.test(as.formula(paste(metric, "~", group_col)), data = df)
  if (kw$p.value < 0.05) {
    dunn <- dunnTest(as.formula(paste(metric, "~", group_col)), data = df, method = "bh")
    cld <- cldList(comparison = dunn$res$Comparison, p.value = dunn$res$P.adj, threshold = 0.05)
    return(cld)
  } else {
    return(NULL)
  }
}

# Get significance letters
letters_observed <- get_signif_letters(alpha_df, "Observed")

# Prepare annotation positions
if (!is.null(letters_observed)) {
  annotation_df <- letters_observed %>%
    rename(Group = Group) %>%
    mutate(y = max(alpha_df$Observed, na.rm = TRUE) * 1.05)
} else {
  annotation_df <- NULL
}

# PLot
p <- ggplot(alpha_df, aes(x = Group, y = Observed, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 3)

# Add letters only if annotation_df exists
if (!is.null(annotation_df)) {
  p <- p + geom_text(data = annotation_df, aes(x = Group, y = y, label = Letter),
                     color = "black", size = 5)
}

p <- p +
  labs(title = "Observed Species Richness", y = "Observed", x = "Group") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank()
  )

print(p)


# Save high-resolution PNG
ggsave("Observed_Richness_AlphaDiversity.png", plot = p, width = 8, height = 6, dpi = 600)


---------------------
# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(FSA)          # For Dunn's Test
library(rcompanion)   # For cldList (compact letter display)
library(ggpubr)
library(patchwork)    # For combining multiple ggplots

# Set working directory
# setwd(".")  # set working directory to the folder containing output.biom and metadata.csv

# Import biom file
merged_metagenomes <- import_biom("output.biom")

# Clean taxonomy labels (remove prefixes like "k__", "p__")
tax_table(merged_metagenomes) <- apply(tax_table(merged_metagenomes), 2, function(x) sub("^[a-z]__", "", x))
colnames(tax_table(merged_metagenomes)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Subset bacteria only
merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria")

# Estimate alpha diversity metrics
alpha_df <- estimate_richness(merged_metagenomes, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
alpha_df$SampleID <- trimws(rownames(alpha_df))  # Trim whitespaces just in case

# Load metadata (do NOT use row.names)
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)
metadata$SampleID <- trimws(metadata$SampleID)  # Ensure no hidden spaces

# Merge alpha diversity with metadata
alpha_df <- merge(alpha_df, metadata, by = "SampleID")

# Function to perform Kruskal-Wallis and post-hoc Dunn test, then assign significance letters
get_signif_letters <- function(df, metric, group_col = "Group") {
  kw <- kruskal.test(as.formula(paste(metric, "~", group_col)), data = df)
  print(paste("Kruskal-Wallis test for", metric))
  print(kw)
  if (kw$p.value < 0.05) {
    dunn <- dunnTest(as.formula(paste(metric, "~", group_col)), data = df, method = "bh")
    cld <- cldList(comparison = dunn$res$Comparison,
                   p.value = dunn$res$P.adj,
                   threshold = 0.05)
    return(cld)
  }
  return(NULL)
}

# Plotting function: boxplot + jitter + significance letters
plot_alpha_diversity_box <- function(df, metric, label) {
  letters_df <- get_signif_letters(df, metric)
  
  if (!is.null(letters_df)) {
    letters_df <- letters_df %>%
      rename(Group = Group) %>%
      mutate(y = max(df[[metric]], na.rm = TRUE) * 1.05)
  }
  
  p <- ggplot(df, aes_string(x = "Group", y = metric, color = "Group", fill = "Group")) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2) +
    labs(title = paste0(label, ". ", metric, " Diversity"), y = metric, x = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, face = "bold", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13, face = "bold")
    )
  
  if (!is.null(letters_df)) {
    p <- p + geom_text(data = letters_df, aes(x = Group, y = y, label = Letter),
                       color = "black", size = 5, inherit.aes = FALSE)
  }
  
  return(p)
}

# Create plots for each metric
p_observed <- plot_alpha_diversity_box(alpha_df, "Observed", "A")
p_chao1    <- plot_alpha_diversity_box(alpha_df, "Chao1", "B")
p_shannon  <- plot_alpha_diversity_box(alpha_df, "Shannon", "C")
p_simpson  <- plot_alpha_diversity_box(alpha_df, "Simpson", "D")

# Combine plots in 2x2 layout using patchwork
combined_plot <- (p_observed | p_chao1) / (p_shannon | p_simpson)

# Show the final plot
print(combined_plot)

# Save high-resolution image
ggsave("Alpha_Diversity_Boxplots_Combined.png", combined_plot, width = 12, height = 10, dpi = 600)


-------
library(dplyr)

# Function to calculate median and range for a given metric
get_median_range <- function(data, metric) {
  data %>%
    group_by(Group) %>%
    summarise(
      Median = round(median(.data[[metric]], na.rm = TRUE), 2),
      Min = round(min(.data[[metric]], na.rm = TRUE), 2),
      Max = round(max(.data[[metric]], na.rm = TRUE), 2)
    ) %>%
    mutate(Metric = metric) %>%
    select(Metric, everything())
}

# List of metrics
metrics <- c("Observed", "Chao1", "Shannon", "Simpson")

# Get summary for all metrics
summary_stats <- lapply(metrics, function(m) get_median_range(alpha_df, m)) %>%
  bind_rows()

print(summary_stats)
