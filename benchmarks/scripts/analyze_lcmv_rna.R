#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(edgeR)
  library(limma)
  library(RUVSeq)
  library(Biobase)
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

# Exploratory cross-study state atlas. Study cannot be added to the model because
# several states are perfectly confounded with study; RUVr mitigates but cannot
# remove that limitation.
meta$state <- factor(meta$broad_condition)
design <- model.matrix(~ state, meta)
y <- calcNormFactors(DGEList(counts = counts), method = "upperquartile")
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmFit(y, design)
residuals_ruv <- residuals(fit, type = "deviance")
set <- newSeqExpressionSet(counts, phenoData = AnnotatedDataFrame(meta))
ruv <- RUVr(set, rownames(counts), k = 1, residuals_ruv)
ruv_meta <- pData(ruv)
if (!"sample" %in% colnames(ruv_meta)) ruv_meta$sample <- rownames(ruv_meta)
write.table(ruv_meta, file.path(out, "ruvr_sample_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ruv_meta$state <- factor(ruv_meta$broad_condition)
dds_ruv <- DESeqDataSetFromMatrix(counts, ruv_meta, design = ~ W_1 + state)
dds_ruv <- DESeq(dds_ruv, quiet = TRUE)
write_matrix(counts(dds_ruv, normalized = TRUE), file.path(out, "ruvr_normalized_counts.tsv"))
broad <- read.delim(file.path(project, "metadata/comparisons_broad.tsv"), stringsAsFactors = FALSE)
broad_dir <- file.path(out, "exploratory_broad_state")
dir.create(broad_dir, recursive = TRUE, showWarnings = FALSE)
for (i in seq_len(nrow(broad))) {
  item <- broad[i, ]
  result <- as.data.frame(results(dds_ruv, contrast = c("state", item$cond1, item$cond2)))
  result <- cbind(gene_key(rownames(result)), result)
  write.table(result, file.path(broad_dir, paste0(item$comparison, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

writeLines(c(
  "Primary RNA differential-expression results use each study's paper-specific k=31 Kallisto layer and within-study DESeq2 contrasts.",
  "Beltra contrasts are also reported using the paper-style TMM/voom/limma route.",
  "The uniform k=21 Kallisto layer is required by the 25-bp Milner reads and supports cross-study visualization and ATAC/RNA correlation.",
  "Guan sensitivity results exclude GSM3045265, whose paper-specific k=31 pseudoalignment rate was 37.4%.",
  "The broad-state RUVr results are exploratory because study and state are partly or perfectly confounded."
), file.path(out, "INTERPRETATION.txt"))
