#Install packages
install.packages("ggpubr")  # if not already installed

# Load required packages
# ========================
library(tidyverse)
library(ggplot2)
library(ggpubr)

# ========================
# Create CFU/mL dataset (4 replicates per group)
# ========================
cfu_data <- data.frame(
  Group = c("CR1", "CR2", "TR1", "TR2", "TR3"),
  Rep1 = c(53, 75, 24, 12, 17),
  Rep2 = c(39, 33, 16, 15, 14),
  Rep3 = c(44, 55, 15, 14, 15),
  Rep4 = c(58, 47, 19, 11, 19)
)

# ========================
# Convert to long format
# ========================
cfu_long <- pivot_longer(cfu_data, cols = starts_with("Rep"),
                         names_to = "Replicate", values_to = "CFU")

# ========================
# Create the plot and assign to variable
# ========================
p <- ggplot(cfu_long, aes(x = Group, y = CFU, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  stat_compare_means(method = "kruskal.test", label.y = 85, label = "p.format") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 15)),
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 15)),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12)
  ) +
  labs(
    y = expression("Average of TVBC × 10"^3 * " (CFU/mL)"),
    x = "Group"
  )

# ========================
# Print the plot
# ========================
print(p)

ggsave("cfu_counts_600dpi.png",
       plot = p, dpi = 600, width = 8, height = 6, units = "in",
       bg = "white")  # add white background

----------------------------------------------------
# Visivble asterik and line
# Create CFU/mL dataset (4 replicates per group)
cfu_data <- data.frame(
  Group = c("CR1", "CR2", "TR1", "TR2", "TR3"),
  Rep1 = c(53, 75, 24, 12, 17),
  Rep2 = c(39, 33, 16, 15, 14),
  Rep3 = c(44, 55, 15, 14, 15),
  Rep4 = c(58, 47, 19, 11, 19)
)

# Convert to long format
cfu_long <- pivot_longer(cfu_data, cols = starts_with("Rep"),
                         names_to = "Replicate", values_to = "CFU")


# Basic bar plot with jittered points
p <- ggplot(cfu_long, aes(x = Group, y = CFU, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "CFU/mL Counts with Statistical Significance",
       y = "CFU/mL")

# Add global Kruskal–Wallis test p-value
p + stat_compare_means(method = "kruskal.test", label.y = 85)  # adjust label.y as needed

p + stat_compare_means(method = "kruskal.test", label.y = 85) +
  stat_compare_means(comparisons = list(c("CR1", "TR1"),
                                        c("CR1", "TR2"),
                                        c("CR2", "TR3")),
                     method = "wilcox.test", label = "p.signif")


ggsave("cfu_significance_plot.png",
       plot = p, width = 6, height = 4, dpi = 600, units = "in")


