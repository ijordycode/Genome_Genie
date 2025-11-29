library(DESeq2)
library(tidyverse)

# Directory containing coutn txt files
count_dir <- "Outputs/counts/"   # <-- change if needed

# Your metadata CSV file
# Must contain columns: sample , condition
sample_csv <- "Outputs/samplesheet.csv"

# -------------------------------
# READ SAMPLE METADATA
# -------------------------------
samples <- read.csv(sample_csv, stringsAsFactors = FALSE)

# Ensure sample column matches EXACT filenames (without .txt)
# Example: sample1.txt -> sample = "sample1"
samples$sample <- as.character(samples$sample)

# Build filenames using the sample column
count_files <- file.path(count_dir, paste0(samples$sample, ".txt"))

# Read all count files
tables <- lapply(count_files, function(x) {
  read.delim(x, comment.char = "#")
})

# Gene IDs from first file
genes <- tables[[1]]$Geneid

# Extract the last column (counts) from each file
count_list <- lapply(tables, function(tbl) tbl[, ncol(tbl)])


countdata <- do.call(cbind, count_list)
rownames(countdata) <- genes
colnames(countdata) <- samples$sample

# -------------------------------
# BUILD DESEQ2 OBJECT
# -------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = samples,
  design = ~ condition
)

# Set rownames of metadata
rownames(colData(dds)) <- colData(dds)$sample

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# -------------------------------
# RUN DESEQ2
# -------------------------------
dds <- DESeq(dds)

# Extract results
res <- results(dds)

write.csv(as.data.frame(res), "Outputs/counts/deseq2_results.csv")