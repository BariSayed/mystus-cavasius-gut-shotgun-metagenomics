# =========================
# Load libraries
# =========================
library(tidyverse)
library(forcats)
library(patchwork)

# =========================
# Data input
# =========================
df <- tribble(
  ~Conditions, ~Bacterial_strains, ~Sensitivity,
  "Control", "Enterobacter cloacae", "Sensitive",
  "Control", "Staphylococcus sciuri", "Sensitive",
  "Control", "Bacillus pumilus", "Sensitive",
  "Control", "Bacillus altitudinis", "Resistant",
  "Control", "Klebsiella pneumoniae", "Sensitive",
  "Control", "Escherichia coli", "Intermediate",
  "Control", "Bacillus tequilensis", "Resistant",
  "Control", "Bacillus atrophaeus", "Intermediate",
  "Control", "Bacillus amyloliquefaciens", "Resistant",
  "Control", "Enterobacter hormaechei", "Sensitive",
  "Treatment", "Bacillus tequilensis", "Resistant",
  "Treatment", "Aeromonas sp.", "Resistant",
  "Treatment", "Aeromonas salmonicida", "Sensitive",
  "Treatment", "Aeromonas caviae", "Resistant",
  "Treatment", "Bacillus aerius", "Resistant",
  "Treatment", "Enterobacter sp.", "Resistant",
  "Treatment", "Enterobacter hormaechei", "Resistant",
  "Treatment", "Bacillus subtilis", "Resistant",
  "Treatment", "Bacillus pumilus", "Resistant",
  "Treatment", "Bacillus cereus", "Sensitive",
  "Treatment", "Bacillus toyonensis", "Sensitive",
  "Treatment", "Enterobacter cloacae", "Resistant",
  "Treatment", "Bacillus thuringiensis", "Resistant",
  "Treatment", "Bacillus sp", "Resistant",
  "Treatment", "Bacillus thuringiensis", "Resistant",
  "Treatment", "Staphylococcus sp.", "Resistant"
)

# =========================
# Summarize data
# =========================
summary_counts <- df %>%
  group_by(Conditions, Sensitivity) %>%
  summarise(Count = n(), .groups = "drop")

summary_counts$Sensitivity <- factor(summary_counts$Sensitivity, levels = c("Sensitive", "Intermediate", "Resistant"))

# =========================
# Create bar plot (Plot A)
# =========================
bar_plot <- ggplot(summary_counts, aes(x = Conditions, y = Count, fill = Sensitivity)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.8), 
            vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Sensitive" = "#1b9e77", 
                               "Intermediate" = "#d95f02", 
                               "Resistant" = "#7570b3")) +
  labs(title = "A. Sensitivity Distribution",
       x = "Conditions",
       y = "Number of Strains",
       fill = "Sensitivity") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )

# =========================
# Prepare italic labels for strains
# =========================
# Create a named vector of expressions for italic scientific names
species_labels <- setNames(
  lapply(unique(df$Bacterial_strains), function(s) bquote(italic(.(s)))),
  unique(df$Bacterial_strains)
)

# =========================
# Create tile plot (Plot B) with italic labels
# =========================
df$Sensitivity <- factor(df$Sensitivity, levels = c("Sensitive", "Intermediate", "Resistant"))

species_tile_plot <- ggplot(df, aes(x = Conditions, y = fct_rev(Bacterial_strains), fill = Sensitivity)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("Sensitive" = "#1b9e77", 
                               "Intermediate" = "#d95f02", 
                               "Resistant" = "#7570b3")) +
  scale_y_discrete(
    labels = species_labels  # Apply italic labels here
  ) +
  labs(title = "B. Species Sensitivity",
       x = "Conditions", y = "Bacterial Strains", fill = "Sensitivity") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),   # Increased & bold species labels
    axis.text.x = element_text(face = "bold", size = 14),   # Increased & bold Control/Treatment labels
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "top"
  )

# =========================
# Combine plots
# =========================
combined_plot <- bar_plot | species_tile_plot

# =========================
# Display combined plot
# =========================
print(combined_plot)

# =========================
# Save combined figure in high resolution
# =========================
ggsave(
  filename = "Combined_Sensitivity_Plots.png",
  plot = combined_plot,
  width = 16,
  height = 10,
  dpi = 600,
  units = "in"
)
