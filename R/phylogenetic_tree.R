# Phylogeny Tree

#Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("treeio", "ggtree"))
install.packages(c("ape", "phangorn", "ggplot2"))  # install only if not yet installed


# Load libraries
library(ape)
library(phangorn)
library(ggtree)
library(treeio)
library(ggplot2)

# Set working directory
# setwd(".")  # set working directory to the folder containing output.biom and metadata.csv

# Read your aligned sequences
alignment <- read.dna("metagenome_36_seq_maft-align.fasta", format = "fasta")

# Convert alignment to phyDat object for ML analysis
phy <- phyDat(alignment, type = "DNA")

# Calculate distance matrix and build initial NJ tree
dm <- dist.ml(phy)
nj_tree <- NJ(dm)

# Optimize maximum likelihood tree using GTR+Gamma model
fit <- pml(nj_tree, data = phy)

# Optimize ML tree with GTR+Gamma model
fit_opt <- optim.pml(fit, model = "GTR", optGamma = TRUE, rearrangement = "stochastic")
ml_tree <- fit_opt$tree

# Perform 1000 bootstrap replicates
bs <- bootstrap.pml(fit_opt, bs=1000, optNni=TRUE, control = pml.control(trace=0))

# Determine maximum x to avoid cutting
max_x <- max(nodeHeights(ml_tree))
cat("Max x-coordinate of the tree:", max_x, "\n")

# Determine maximum x to avoid cutting tips
tip_depths <- ape::node.depth.edgelength(ml_tree)
max_x <- max(tip_depths)
cat("Max x-coordinate of the tree:", max_x, "\n")


# Create a bootstrap-annotated tree plot, then add better styling and squeeze
p <- ggtree(plotBS(ml_tree, bs, p=50)) +   # adds bootstrap values >=50%
  geom_tiplab(size=7) +                    # bigger tip labels
  theme_tree2() +
  geom_treescale(x=0, y=0, width=0.02) +
  xlim_tree(1)                            # squeeze tree horizontally into x=0–1

print(p)

# Save the tree plot to a PNG file
ggsave("ML_phylogenetic_tree_bootstrap.svg", plot = p, width=20, height=14, dpi=600)
ggsave("ML_phylogenetic_tree_bootstrap.svg", plot = p, width=20, height=14, dpi=600)
ggsave("ML_phylogenetic_tree_bootstrap.pdf", plot = p, width=20, height=14, dpi=600)

# Save the tree itself to a Newick file in your desired directory
write.tree(ml_tree, file = "ML_phylogenetic_tree.nwk")
