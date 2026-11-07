
# ====== Load Libraries ======
library(phyloseq)
library(dplyr)
library(tidyr)
library(openxlsx)

# ====== Import Data ======
physeq_raw <- import_biom("output.biom")
metadata <- read.csv("metadata.csv", row.names = 1, check.names = FALSE)
sample_data(physeq_raw) <- sample_data(metadata)

# ====== Clean Taxonomy ======
tax_table(physeq_raw) <- apply(tax_table(physeq_raw), 2, function(x) sub("^[a-z]__", "", x))
colnames(tax_table(physeq_raw)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ====== Filter Bacteria (excluding Homo sapiens) ======
physeq_filtered <- subset_taxa(physeq_raw, Kingdom == "Bacteria" & !grepl("Homo sapiens", Species, ignore.case = TRUE))

# ====== Calculate relative abundances ======
physeq_rel <- transform_sample_counts(physeq_filtered, function(x) x / sum(x))

# ====== Extract raw counts table ======
otu_count_mat <- as(otu_table(physeq_filtered), "matrix") %>% as.data.frame()
otu_count_mat$OTU_ID <- rownames(otu_count_mat)

# ====== Extract relative abundance table ======
otu_rel_mat <- as(otu_table(physeq_rel), "matrix") %>% as.data.frame()
otu_rel_mat$OTU_ID <- rownames(otu_rel_mat)

# ====== Extract taxonomy table ======
tax_df <- as(tax_table(physeq_filtered), "matrix") %>% as.data.frame()
tax_df$OTU_ID <- rownames(tax_df)

# ====== Merge raw counts with taxonomy ======
otu_count_tax <- merge(otu_count_mat, tax_df, by = "OTU_ID", all.x = TRUE)
otu_count_tax$Species_Name <- paste(otu_count_tax$Genus, otu_count_tax$Species, sep = " ")

# ====== Merge relative abundances with taxonomy ======
otu_rel_tax <- merge(otu_rel_mat, tax_df, by = "OTU_ID", all.x = TRUE)
otu_rel_tax$Species_Name <- paste(otu_rel_tax$Genus, otu_rel_tax$Species, sep = " ")

# ====== View actual columns/sample names in counts table ======
cat("Available columns (samples) in counts table:\n")
print(colnames(otu_count_tax))

# ====== Select desired columns with correct sample names ======
desired_samples <- c("Control_S1", "Control_S8", "SingleAbx_S2", "SingleAbx_S3", "SingleAbx_S4")
desired_cols <- c("OTU_ID", "Species_Name", desired_samples)

existing_cols_count <- desired_cols[desired_cols %in% colnames(otu_count_tax)]
existing_cols_rel <- desired_cols[desired_cols %in% colnames(otu_rel_tax)]

otu_count_clean <- otu_count_tax[, existing_cols_count]
otu_rel_clean <- otu_rel_tax[, existing_cols_rel]

# ====== Write to Excel with multiple sheets ======
output_excel <- "otu_table_counts_and_relative.xlsx"
wb <- createWorkbook()
addWorksheet(wb, "Raw_Counts")
writeData(wb, "Raw_Counts", otu_count_clean)
addWorksheet(wb, "Relative_Abundance")
writeData(wb, "Relative_Abundance", otu_rel_clean)
saveWorkbook(wb, output_excel, overwrite = TRUE)

cat(paste0("✅ Excel file '", output_excel, "' created successfully with raw counts and relative abundances.\n"))





---------------------------
# Load necessary packages
library(phyloseq)
library(dplyr)

# Assume your phyloseq object is physeq_filtered (already loaded and filtered)
# If not, replace physeq_filtered with your object name

# Extract taxonomy table as a data frame
tax_df <- as(tax_table(physeq_filtered), "matrix") %>% as.data.frame(stringsAsFactors = FALSE)

# Define function to count unique non-missing taxa
count_unique <- function(x) {
  length(unique(x[!is.na(x) & x != ""]))
}

# Count unique taxa at each rank
num_phyla <- count_unique(tax_df$Phylum)
num_classes <- count_unique(tax_df$Class)
num_families <- count_unique(tax_df$Family)
num_genera <- count_unique(tax_df$Genus)
num_species <- count_unique(tax_df$Species)

# Print formatted summary sentence
cat(sprintf("A total of %d phyla, %d classes, %d families, %d genera, and %d bacterial taxa at the species level were identified in the gut of Mystus cavasius.\n",
            num_phyla, num_classes, num_families, num_genera, num_species))
