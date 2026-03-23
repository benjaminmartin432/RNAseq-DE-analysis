# ============================================================
# Transcriptomic Analysis of Opaque-2 Loss-of-Function
# in Zea mays
#
# Author: Ben Martin
# Data:   EMBL-EBI Expression Atlas, Accession E-CURD-41
#         Zhan et al. (2018) Opaque-2 Regulates a Complex
#         Gene Network Associated with Cell Differentiation
#         and Storage Functions of Maize Endosperm.
# ===========================================================

# ============================================================
# SECTION 1: LOAD OR INSTALL LIBRARIES
# ============================================================

requiredPackages = c('BiocManager', 'dplyr', 'stringr', 'tidyr')
for (p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)}
BiocManagerPKG <- c('DESeq2', 'pheatmap', 'EnhancedVolcano')
for (p in BiocManagerPKG){
  if(!require(p,character.only = TRUE)) BiocManager::install(p)
  library(p, character.only = TRUE)}

# ============================================================
# USER CONFIGURATION
# Set your preferences here before running the script.
# ============================================================

# Set to TRUE to save plots and results to disk, FALSE to skip saving
SAVE_OUTPUT <- FALSE

# If SAVE_OUTPUT is TRUE, set your desired output directories here.
# The directories will be created automatically if they do not exist.
# Use "." to save to the current working directory.
PLOT_DIR    <- "plots"
RESULTS_DIR <- "results"

# Create output directories if saving is enabled
if (SAVE_OUTPUT) {
  dir.create(PLOT_DIR,    showWarnings = FALSE, recursive = TRUE)
  dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
}

# Helper functions: return full path if saving, NULL otherwise
plot_path    <- function(filename) if (SAVE_OUTPUT) file.path(PLOT_DIR,    filename) else NULL
results_path <- function(filename) if (SAVE_OUTPUT) file.path(RESULTS_DIR, filename) else NULL

# ============================================================
# SECTION 2: LOAD DATA
# ============================================================

# These files are all present in the repository within the data folder and can be downloaded directly if desired 

# --- Raw count data ---
# --- If needed this file can be downloaded directly from the repository under E-CURD-41-raw-counts.tsv ---
# 44,303 genes examined across 6 samples (3 wildtype, 3 mutant)
counts <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-CURD-41/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")

# --- Sample metadata ---
# --- if needed this file can be downloaded directly from the repository under E-CURD-41-experiment-design.tsv ---
metadata <- read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-CURD-41/resources/ExperimentDesignFile.RnaSeq/experiment-design")

# --- Gene annotation ---
# --- If needed this file can be downloaded directly from the repository under Zm00001eb.1.fulldata.txt ---
# Downloaded from MaizeGDB (Zm-B73-REFERENCE-NAM-5.0)
download.file(
  "https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm00001eb.1.fulldata.txt.gz",
  destfile = "maize.txt.gz"
)

genes <- read.delim(
  gzfile("maize.txt.gz"),
  comment.char = "#",
  header        = FALSE,
  sep           = "\t"
)



# ============================================================
# SECTION 3: GENE NAME MAPPING
# Maps Zm00001eb IDs to readable gene names where available,
# falling back to alternative IDs or the gene ID itself
# ============================================================

gene_map <- data.frame(
  gene_id            = genes$V2,
  gene_name_raw      = genes$V12,
  gene_name_fallback = genes$V11,
  stringsAsFactors   = FALSE
)

# Encode placeholder dashes as NA
gene_map$gene_name_raw[gene_map$gene_name_raw == "-"]           <- NA
gene_map$gene_name_fallback[gene_map$gene_name_fallback == "-"] <- NA

# Priority: full name > alternative ID > gene ID
gene_map$gene_name <- ifelse(
  !is.na(gene_map$gene_name_raw),
  gene_map$gene_name_raw,
  ifelse(
    !is.na(gene_map$gene_name_fallback),
    gene_map$gene_name_fallback,
    gene_map$gene_id
  )
)

# Retain only gene ID and resolved name; drop duplicate IDs
gene_map_final <- gene_map[, c("gene_id", "gene_name")]
gene_map_final <- gene_map_final[!duplicated(gene_map_final$gene_id), ]


# ============================================================
# SECTION 4: FORMAT DATA FOR DESEQ2
# ============================================================

# DESeq2 expects gene IDs as count matrix rownames
rownames(counts) <- counts$Gene.ID
counts           <- counts[-c(1, 2)]  # remove Gene.ID and Gene.Name columns

# DESeq2 expects sample IDs as metadata rownames
rownames(metadata) <- metadata$Run
metadata           <- metadata[colnames(counts), ]

# Retain only genotype column
metadata <- metadata[, "Factor.Value.genotype.", drop = FALSE]
colnames(metadata) <- "genotype"

# Recode genotype labels to clean factor levels
metadata$genotype[metadata$genotype == "wild type genotype"]       <- "wildtype"
metadata$genotype[metadata$genotype == "O2 loss of function mutant"] <- "mutant"
metadata$genotype <- factor(metadata$genotype, levels = c("wildtype", "mutant"))


# ============================================================
# SECTION 5: DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS
# Design: ~ genotype (wildtype as reference level)
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ genotype
)

# Remove genes with very low total counts across all samples
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res    <- results(dds)
res.df <- as.data.frame(res)

# Add gene ID as a column (rownames are gene IDs)
res.df$gene_id <- rownames(res.df)

# Map gene names into results dataframe
res.df$gene_name <- gene_map_final$gene_name[
  match(res.df$gene_id, gene_map_final$gene_id)
]

cat("Total genes tested:", nrow(res.df), "\n")
cat("Significant DE genes (padj < 0.05):",
    sum(!is.na(res.df$padj) & res.df$padj < 0.05), "\n")


# ============================================================
# SECTION 6: VARIANCE STABILISING TRANSFORMATION
# Used for exploratory visualisation (PCA and heatmap).
# blind = FALSE uses the experimental design to inform
# the transformation, appropriate for post-QC visualisation.
# ============================================================

vsd <- vst(dds, blind = FALSE)


# ============================================================
# SECTION 7: PCA
# Assesses replicate clustering and confirms genotype as
# the dominant source of transcriptional variance
# ============================================================

pcaData   <- plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = genotype)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

print(pca_plot)

if (SAVE_OUTPUT) {
  ggsave(plot_path("PCA_plot.png"),
         plot   = pca_plot,
         width  = 7,
         height = 6,
         units  = "in",
         dpi    = 300)
}


# ============================================================
# SECTION 8: MA PLOT
# Log-ratio vs mean expression across all genes
# ============================================================

# Always print to console
plotMA(res)

# Optionally save to file
if (SAVE_OUTPUT) {
  png(plot_path("MA_plot.png"), width = 8, height = 6, units = "in", res = 300)
  plotMA(res)
  dev.off()
}


# ============================================================
# SECTION 9: VOLCANO PLOT
# log2FoldChange vs -log10(pvalue), labelled with gene names
# ============================================================

volcano_plot <- EnhancedVolcano(
  res.df,
  lab = res.df$gene_name,
  x   = "log2FoldChange",
  y   = "pvalue"
)

print(volcano_plot)

if (SAVE_OUTPUT) {
  ggsave(plot_path("volcano_plot.png"),
         plot   = volcano_plot,
         width  = 12,
         height = 10,
         units  = "in",
         dpi    = 300)
}


# ============================================================
# SECTION 10: HEATMAP — TOP 20 DE GENES
# Row-scaled VST expression values, hierarchical clustering
# of both genes and samples
# ============================================================

# Select the 20 genes with the smallest adjusted p-values
top20 <- head(order(res.df$padj, na.last = TRUE), 20)
mat   <- assay(vsd)[rownames(res.df)[top20], ]

# Substitute gene names for row labels
name_map     <- setNames(res.df$gene_name, res.df$gene_id)
rownames(mat) <- name_map[match(rownames(mat), names(name_map))]

# Row-scale so colour reflects relative change, not absolute expression
mat_scaled <- t(scale(t(mat)))

# Always print to console (filename argument omitted)
pheatmap(mat_scaled,
         cluster_rows   = TRUE,
         cluster_cols   = TRUE,
         show_rownames  = TRUE,
         annotation_col = metadata,
         fontsize_row   = 10,
         fontsize_col   = 10,
         color          = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
         breaks         = seq(-2, 2, length.out = 101),
         main           = "Top 20 DE Genes: Opaque-2 vs Wildtype",
         border_color   = NA,
         treeheight_row = 40,
         treeheight_col = 20)

# Optionally save to file — pheatmap suppresses console display
# when filename is set, so this call is kept separate
if (SAVE_OUTPUT) {
  pheatmap(mat_scaled,
           cluster_rows   = TRUE,
           cluster_cols   = TRUE,
           show_rownames  = TRUE,
           annotation_col = metadata,
           fontsize_row   = 10,
           fontsize_col   = 10,
           color          = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
           breaks         = seq(-2, 2, length.out = 101),
           main           = "Top 20 DE Genes: Opaque-2 vs Wildtype",
           border_color   = NA,
           treeheight_row = 40,
           treeheight_col = 20,
           filename       = plot_path("heatmap_top20.png"),
           width          = 8,
           height         = 10)
}


# ============================================================
# SECTION 11: GO TERM PARSING
# Extracts GO:XXXXXXX IDs and descriptions from the MaizeGDB
# annotation. Splits on the GO: prefix to avoid the
# embedded-comma problem within GO term descriptions.
# ============================================================

gene_go_raw <- genes[, c("V2", "V14")]
colnames(gene_go_raw) <- c("gene_id", "go_terms")
gene_go_raw <- gene_go_raw[
  !is.na(gene_go_raw$go_terms) & gene_go_raw$go_terms != "", 
]

go_long <- lapply(seq_len(nrow(gene_go_raw)), function(i) {
  gid   <- gene_go_raw$gene_id[i]
  terms <- gene_go_raw$go_terms[i]
  
  go_ids   <- stringr::str_extract_all(terms, "GO:\\d{7}")[[1]]
  go_descs <- stringr::str_extract_all(
    terms, "(?<=GO:\\d{7}=)[^|]*?(?=,?GO:\\d{7}=|$)"
  )[[1]]
  go_descs <- stringr::str_remove(go_descs, ",$")
  go_descs <- stringr::str_trim(go_descs)
  
  n <- min(length(go_ids), length(go_descs))
  if (n == 0) return(NULL)
  
  data.frame(
    gene_id          = gid,
    go_id            = go_ids[1:n],
    go_desc          = go_descs[1:n],
    stringsAsFactors = FALSE
  )
})

go_long <- do.call(rbind, go_long)
go_long <- go_long[!duplicated(go_long[, c("gene_id", "go_id")]), ]

cat("Genes with GO annotations:", length(unique(go_long$gene_id)), "\n")
cat("Unique GO terms found:",     length(unique(go_long$go_id)),   "\n")


# ============================================================
# SECTION 12: IDENTIFY SIGNIFICANTLY DE GENES
# ============================================================

sig_genes             <- res.df[!is.na(res.df$padj) & res.df$padj < 0.05, ]
sig_genes$direction   <- ifelse(
  sig_genes$log2FoldChange > 0, "upregulated", "downregulated"
)

cat("Upregulated:", sum(sig_genes$direction == "upregulated"), "\n")
cat("Downregulated:", sum(sig_genes$direction == "downregulated"), "\n")

# Join DE genes with GO terms
sig_go <- merge(
  sig_genes[, c("gene_id", "log2FoldChange", "padj", "direction")],
  go_long,
  by = "gene_id"
)

cat("Significant genes with GO annotations:",
    length(unique(sig_go$gene_id)), "\n")


# ============================================================
# SECTION 13: GO TERM FREQUENCY TABLE
# ============================================================

go_freq <- sig_go |>
  dplyr::group_by(go_id, go_desc) |>
  dplyr::summarise(
    gene_count = dplyr::n_distinct(gene_id),
    up_count   = dplyr::n_distinct(gene_id[direction == "upregulated"]),
    down_count = dplyr::n_distinct(gene_id[direction == "downregulated"]),
    .groups    = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(gene_count))

if (SAVE_OUTPUT) {
  write.csv(go_freq, results_path("GO_term_frequency_table.csv"), row.names = FALSE)
}

# ============================================================
# SECTION 14: OVER-REPRESENTATION ANALYSIS (ORA)
# One-sided Fisher's exact test per GO term, comparing
# frequency among DE genes vs whole-genome background.
# P-values corrected with Benjamini-Hochberg method.
# ============================================================

all_annotated_genes <- length(unique(go_long$gene_id))
total_sig           <- length(unique(sig_go$gene_id))

go_background <- go_long |>
  dplyr::group_by(go_id, go_desc) |>
  dplyr::summarise(
    background_count = dplyr::n_distinct(gene_id),
    .groups          = "drop"
  )

ora_input <- merge(go_freq, go_background, by = c("go_id", "go_desc"))

fisher_results <- lapply(seq_len(nrow(ora_input)), function(i) {
  a <- ora_input$gene_count[i]
  b <- ora_input$background_count[i] - a
  c <- total_sig - a
  d <- all_annotated_genes - ora_input$background_count[i] - c
  
  if (any(c(a, b, c, d) < 0)) {
    return(data.frame(p_value = NA, odds_ratio = NA))
  }
  
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  data.frame(p_value = ft$p.value, odds_ratio = as.numeric(ft$estimate))
})

fisher_df          <- do.call(rbind, fisher_results)
ora_results        <- cbind(ora_input, fisher_df)
ora_results$p_adj  <- p.adjust(ora_results$p_value, method = "BH")
ora_results        <- ora_results[order(ora_results$p_adj), ]

cat("\nSignificantly enriched GO terms (padj < 0.05):",
    sum(!is.na(ora_results$p_adj) & ora_results$p_adj < 0.05), "\n")

if (SAVE_OUTPUT) {
  write.csv(ora_results, results_path("GO_ORA_results.csv"), row.names = FALSE)
}


# ============================================================
# SECTION 15: GO TERM VISUALISATIONS
# ============================================================

top_n <- 20

# --- Plot 1: Stacked bar chart — GO term frequency ---
go_freq_top         <- head(go_freq, top_n)
go_freq_top$go_label <- stringr::str_wrap(
  paste0(go_freq_top$go_id, ": ", go_freq_top$go_desc), width = 45
)

bar_data <- rbind(
  data.frame(go_label   = go_freq_top$go_label,
             gene_count = go_freq_top$gene_count,
             count      = go_freq_top$up_count,
             direction  = "Upregulated"),
  data.frame(go_label   = go_freq_top$go_label,
             gene_count = go_freq_top$gene_count,
             count      = go_freq_top$down_count,
             direction  = "Downregulated")
)

p1 <- ggplot(bar_data,
             aes(x = reorder(go_label, gene_count),
                 y = count,
                 fill = direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c("Upregulated" = "#d73027", "Downregulated" = "#4575b4")
  ) +
  coord_flip() +
  labs(
    title    = "Top GO Terms in Differentially Expressed Genes",
    subtitle = "Opaque-2 Loss-of-Function vs Wildtype",
    x        = NULL,
    y        = "Number of DE Genes",
    fill     = "Direction"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p1)

if (SAVE_OUTPUT) {
  ggsave(plot_path("GO_barplot_frequency.png"),
         plot   = p1,
         width  = 12,
         height = 8,
         units  = "in",
         dpi    = 300)
}

# --- Plot 2: ORA dot plot — odds ratio and significance ---
ora_sig            <- ora_results[!is.na(ora_results$p_adj) & ora_results$p_adj < 0.05, ]
ora_sig            <- head(ora_sig[order(-ora_sig$odds_ratio), ], top_n)
ora_sig$go_label   <- stringr::str_wrap(
  paste0(ora_sig$go_id, ": ", ora_sig$go_desc), width = 50
)
ora_sig$log10_padj <- -log10(ora_sig$p_adj)

p2 <- ggplot(ora_sig,
             aes(x      = odds_ratio,
                 y      = reorder(go_label, odds_ratio),
                 size   = gene_count,
                 colour = log10_padj)) +
  geom_point() +
  scale_colour_gradient(
    low  = "#fee090",
    high = "#d73027",
    name = "-log10(adj. p)"
  ) +
  scale_size_continuous(name = "DE Gene Count", range = c(3, 10)) +
  labs(
    title    = "GO Term Over-Representation Analysis",
    subtitle = "Fisher's Exact Test, Benjamini-Hochberg correction",
    x        = "Odds Ratio",
    y        = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p2)

if (SAVE_OUTPUT) {
  ggsave(plot_path("GO_dotplot_ORA.png"),
         plot   = p2,
         width  = 13,
         height = 8,
         units  = "in",
         dpi    = 300)
}

# --- Plot 3: Directional dot plot — up vs downregulated ---
up_go           <- sig_go[sig_go$direction == "upregulated", ]
up_go           <- aggregate(gene_id ~ go_id + go_desc, data = up_go,
                             FUN = function(x) length(unique(x)))
colnames(up_go)[3] <- "gene_count"
up_go$direction <- "Upregulated"

dn_go           <- sig_go[sig_go$direction == "downregulated", ]
dn_go           <- aggregate(gene_id ~ go_id + go_desc, data = dn_go,
                             FUN = function(x) length(unique(x)))
colnames(dn_go)[3] <- "gene_count"
dn_go$direction <- "Downregulated"

dir_go <- rbind(up_go, dn_go)
dir_go <- merge(dir_go,
                ora_results[, c("go_id", "p_adj", "odds_ratio")],
                by = "go_id")
dir_go <- dir_go[!is.na(dir_go$p_adj) & dir_go$p_adj < 0.05, ]

dir_go <- do.call(rbind, lapply(split(dir_go, dir_go$direction), function(x) {
  head(x[order(-x$gene_count), ], 15)
}))

dir_go$go_label <- stringr::str_wrap(
  paste0(dir_go$go_id, ": ", dir_go$go_desc), width = 45
)

p3 <- ggplot(dir_go,
             aes(x      = direction,
                 y      = reorder(go_label, gene_count),
                 size   = gene_count,
                 colour = -log10(p_adj))) +
  geom_point() +
  scale_colour_gradient(
    low  = "#74add1",
    high = "#d73027",
    name = "-log10(adj. p)"
  ) +
  scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
  labs(
    title    = "Enriched GO Terms by Direction of Regulation",
    subtitle = "Opaque-2 Loss-of-Function vs Wildtype",
    x        = NULL,
    y        = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p3)

if (SAVE_OUTPUT) {
  ggsave(plot_path("GO_dotplot_directional.png"),
         plot   = p3,
         width  = 13,
         height = 9,
         units  = "in",
         dpi    = 300)
}