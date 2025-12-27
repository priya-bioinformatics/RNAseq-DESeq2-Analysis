# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c(
  "ExperimentHub",
  "airway",
  "DESeq2",
  "SummarizedExperiment",
  "EnhancedVolcano",
  "pheatmap",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot"
))

# Load libraries & data
library(ExperimentHub)
library(airway)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)

# Load airway dataset
data("airway")

# Prepare DESeq2 dataset
dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex   # control for cell, test dex treatment
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Remove NA padj values
res <- res[!is.na(res$padj), ]

# Order by adjusted p-value
res <- res[order(res$padj), ]

# View top genes
head(res)

# Save differential expression results
write.csv(
  as.data.frame(res),
  file = "airway_DESeq2_results.csv"
)

# MA plot
plotMA(res,
       main = "DESeq2 MA Plot – Airway",
       ylim = c(-5, 5))

# Dispersion plot
plotDispEsts(dds)

# Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = FALSE)

# PCA plot
plotPCA(vsd, intgroup = "dex")

# Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0,
                title = "Volcano Plot – Airway (Dex vs Control)"
)

# Heatmap of top 20 significant genes
top_genes <- head(rownames(res), 20)

# Extract VST-normalized counts
mat <- assay(vsd)[top_genes, ]

# Column annotation
annotation_col <- as.data.frame(colData(airway)[, "dex", drop = FALSE])
colnames(annotation_col) <- "Treatment"

pheatmap(
  mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  annotation_col = annotation_col,
  main = "Top 20 Differentially Expressed Genes"
)

# GO & KEGG PATHWAY ENRICHMENT

# Select Significant genes
sig_genes <- rownames(res)[res$padj < 0.05]
length(sig_genes)

# Convert Gene Symbols → Entrez IDs
gene_ids <- bitr(
  sig_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# REMOVE duplicated Entrez IDs
gene_ids_unique <- gene_ids[!duplicated(gene_ids$ENTREZID), ]

# GO Enrichment
go_enrich <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Dot plot
dotplot(go_enrich, showCategory = 15) +
  ggtitle("GO Biological Processes – Dexamethasone Response")

# KEGG Pathway Enrichment
kegg_enrich <- enrichKEGG(
  gene         = gene_ids$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Make gene IDs readable
kegg_enrich <- setReadable(
  kegg_enrich,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

# Dot plot
dotplot(kegg_enrich, showCategory = 10) +
  ggtitle("KEGG Pathways – Steroid / Inflammatory Response")

# Save Enrichment Results
write.csv(as.data.frame(go_enrich),
          "GO_enrichment_results.csv")

write.csv(as.data.frame(kegg_enrich),
          "KEGG_enrichment_results.csv")

# Final Checks

# Size factors (normalization)
sizeFactors(dds)

# Sample metadata
colData(dds)
