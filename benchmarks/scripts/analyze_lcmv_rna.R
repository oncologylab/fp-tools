#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(edgeR)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: analyze_lcmv_rna.R PROJECT")
project <- normalizePath(args[[1]], mustWork = TRUE)
metadata_file <- file.path(project, "rna/metadata/samples.tsv")
tx2gene_file <- file.path(project, "rna/metadata/tx2gene.tsv")
out <- file.path(project, "rna/deseq2")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
txmap <- read.delim(tx2gene_file, stringsAsFactors = FALSE)
tx2gene <- unique(txmap[, c("transcript_id", "gene_id")])
symbol_map <- setNames(txmap$gene_symbol, txmap$gene_id)

import_kallisto <- function(samples, layer) {
  files <- file.path(project, "rna", layer, samples, "abundance.tsv")
  names(files) <- samples
  if (any(!file.exists(files))) stop("Missing Kallisto output: ", paste(files[!file.exists(files)], collapse = ", "))
  tximport(
    files,
    type = "kallisto",
    tx2gene = tx2gene,
    dropInfReps = TRUE,
    ignoreAfterBar = TRUE,
    countsFromAbundance = "lengthScaledTPM"
  )
}

gene_key <- function(ids) {
  data.frame(gene_id = ids, gene_symbol = unname(symbol_map[ids]), row.names = ids, check.names = FALSE)
}

write_matrix <- function(matrix, path) {
  key <- gene_key(rownames(matrix))
  result <- cbind(key, as.data.frame(matrix, check.names = FALSE))
  write.table(result, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Uniform k=21 layer: the shorter k-mer is required for the 25-bp Milner reads.
uniform_txi <- import_kallisto(meta$sample, "uniform_kallisto")
uniform_counts <- round(uniform_txi$counts)
storage.mode(uniform_counts) <- "integer"
uniform_counts_dir <- file.path(project, "rna/counts/uniform_kallisto")
dir.create(uniform_counts_dir, recursive = TRUE, showWarnings = FALSE)
write_matrix(uniform_counts, file.path(uniform_counts_dir, "gene_counts_tximport_length_scaled.tsv"))

keep <- rowSums(uniform_counts >= 10) >= 2
counts <- uniform_counts[keep, , drop = FALSE]
meta$condition <- factor(meta$condition)
dds_all <- DESeqDataSetFromMatrix(counts, meta, design = ~ condition)
dds_all <- estimateSizeFactors(dds_all)
write_matrix(counts(dds_all, normalized = TRUE), file.path(out, "uniform_normalized_counts.tsv"))
vst_all <- vst(dds_all, blind = TRUE)
write_matrix(assay(vst_all), file.path(out, "uniform_vst.tsv"))
write.table(
  cor(assay(vst_all), method = "pearson"),
  file.path(out, "descriptive_sample_correlations.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)
pca <- prcomp(t(assay(vst_all)))
pca_rows <- cbind(meta[rownames(pca$x), c("sample", "author", "condition", "broad_condition", "collection")], pca$x[, 1:5, drop = FALSE])
write.table(pca_rows, file.path(out, "descriptive_pca.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Paper-specific k=31 Kallisto layer for the Guan and Beltra studies.
paper_counts <- list()
paper_count_dir <- file.path(project, "rna/counts/paper_specific_tximport")
dir.create(paper_count_dir, recursive = TRUE, showWarnings = FALSE)
for (author in c("Guan", "Beltra")) {
  samples <- meta$sample[meta$author == author]
  imported <- import_kallisto(samples, "paper_specific")
  matrix <- round(imported$counts)
  storage.mode(matrix) <- "integer"
  paper_counts[[author]] <- matrix
  write_matrix(matrix, file.path(paper_count_dir, paste0(tolower(author), "_gene_counts_tximport_length_scaled.tsv")))
}

# Paper-specific TopHat2/HTSeq counts for the two short-read studies.
for (author in c("Milner", "Scott-Browne")) {
  key <- tolower(gsub("-", "_", author))
  path <- file.path(project, "rna/counts/paper_specific", paste0(key, "_gene_counts.tsv"))
  table <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(table) <- table$gene_id
  samples <- meta$sample[meta$author == author]
  missing <- setdiff(samples, colnames(table))
  if (length(missing)) stop("Missing paper-specific HTSeq samples: ", paste(missing, collapse = ", "))
  matrix <- round(as.matrix(table[, samples, drop = FALSE]))
  storage.mode(matrix) <- "integer"
  paper_counts[[author]] <- matrix
}

# Primary inference: paper-specific counts and only within-study contrasts.
comparisons <- read.delim(file.path(project, "metadata/rna_comparisons.tsv"), stringsAsFactors = FALSE)
fine_dir <- file.path(out, "within_study")
dir.create(fine_dir, recursive = TRUE, showWarnings = FALSE)
for (i in seq_len(nrow(comparisons))) {
  item <- comparisons[i, ]
  selected <- meta$condition %in% c(item$cond1, item$cond2)
  submeta <- droplevels(meta[selected, , drop = FALSE])
  author <- unique(submeta$author)
  if (length(author) != 1 || !author %in% names(paper_counts)) stop("Invalid within-study comparison: ", item$comparison)
  subcounts <- paper_counts[[author]][, submeta$sample, drop = FALSE]
  subkeep <- rowSums(subcounts >= 10) >= 2
  dds <- DESeqDataSetFromMatrix(subcounts[subkeep, , drop = FALSE], submeta, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  result <- as.data.frame(results(dds, contrast = c("condition", item$cond1, item$cond2)))
  result <- cbind(gene_key(rownames(result)), result)
  write.table(result, file.path(fine_dir, paste0(item$comparison, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

# Supporting RNA-only contrasts remain separate from paired inference.
supporting <- read.delim(file.path(project, "metadata/supporting_rna_comparisons.tsv"), stringsAsFactors = FALSE)
supporting_dir <- file.path(out, "supporting_rna_only")
dir.create(supporting_dir, recursive = TRUE, showWarnings = FALSE)
for (i in seq_len(nrow(supporting))) {
  item <- supporting[i, ]
  selected <- meta$condition %in% c(item$cond1, item$cond2)
  submeta <- droplevels(meta[selected, , drop = FALSE])
  author <- unique(submeta$author)
  if (length(author) != 1 || !author %in% names(paper_counts)) stop("Invalid supporting comparison: ", item$comparison)
  subcounts <- paper_counts[[author]][, submeta$sample, drop = FALSE]
  subkeep <- rowSums(subcounts >= 10) >= 2
  dds <- DESeqDataSetFromMatrix(subcounts[subkeep, , drop = FALSE], submeta, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  result <- as.data.frame(results(dds, contrast = c("condition", item$cond1, item$cond2)))
  result <- cbind(gene_key(rownames(result)), result)
  write.table(result, file.path(supporting_dir, paste0(item$comparison, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

# Guan sensitivity analysis excluding the naïve library with <50% k=31
# pseudoalignment (GSM3045265). The other two naïve replicates are retained.
guan_sensitivity_dir <- file.path(out, "sensitivity_exclude_GSM3045265")
dir.create(guan_sensitivity_dir, recursive = TRUE, showWarnings = FALSE)
for (i in which(startsWith(comparisons$comparison, "guan_"))) {
  item <- comparisons[i, ]
  selected <- meta$condition %in% c(item$cond1, item$cond2) & meta$sample != "GSM3045265"
  submeta <- droplevels(meta[selected, , drop = FALSE])
  subcounts <- paper_counts[["Guan"]][, submeta$sample, drop = FALSE]
  subkeep <- rowSums(subcounts >= 10) >= 2
  dds <- DESeqDataSetFromMatrix(subcounts[subkeep, , drop = FALSE], submeta, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  result <- as.data.frame(results(dds, contrast = c("condition", item$cond1, item$cond2)))
  result <- cbind(gene_key(rownames(result)), result)
  write.table(result, file.path(guan_sensitivity_dir, paste0(item$comparison, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

# Beltra paper-style TMM/voom/limma sensitivity analysis.
beltra_meta <- droplevels(meta[meta$author == "Beltra", , drop = FALSE])
beltra_counts <- paper_counts[["Beltra"]][, beltra_meta$sample, drop = FALSE]
beltra_group <- factor(beltra_meta$condition)
y_beltra <- DGEList(counts = beltra_counts)
y_beltra <- y_beltra[filterByExpr(y_beltra, group = beltra_group), , keep.lib.sizes = FALSE]
y_beltra <- calcNormFactors(y_beltra, method = "TMM")
design_beltra <- model.matrix(~ 0 + beltra_group)
colnames(design_beltra) <- levels(beltra_group)
voom_beltra <- voom(y_beltra, design_beltra, plot = FALSE)
fit_beltra <- lmFit(voom_beltra, design_beltra)
limma_dir <- file.path(out, "beltra_tmm_voom_limma")
dir.create(limma_dir, recursive = TRUE, showWarnings = FALSE)
for (i in which(startsWith(comparisons$comparison, "beltra_"))) {
  item <- comparisons[i, ]
  contrast <- makeContrasts(contrasts = paste0(item$cond1, "-", item$cond2), levels = design_beltra)
  fit2 <- eBayes(contrasts.fit(fit_beltra, contrast), robust = TRUE)
  result <- topTable(fit2, number = Inf, sort.by = "none")
  result <- cbind(gene_key(rownames(result)), result)
  write.table(result, file.path(limma_dir, paste0(item$comparison, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

writeLines(c(
  "Primary RNA differential-expression results use paper-specific within-study count layers and DESeq2 contrasts.",
  "Guan and Beltra use k=31 Kallisto with tximport lengthScaledTPM counts; Milner and Scott-Browne use TopHat2/HTSeq gene counts.",
  "Beltra contrasts are also reported using the paper-style TMM/voom/limma route.",
  "The uniform k=21 Kallisto layer is required by the 25-bp Milner reads and supports descriptive PCA, correlations, and ATAC/RNA integration only.",
  "Guan sensitivity results exclude GSM3045265, whose paper-specific k=31 pseudoalignment rate was 37.4%.",
  "No pooled cross-study differential test is produced because study and state are partly or perfectly confounded.",
  "Supporting RNA-only contrasts are written separately and are not paired ATAC/RNA evidence."
), file.path(out, "INTERPRETATION.txt"))
