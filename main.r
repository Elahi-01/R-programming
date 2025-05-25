# ==============================
# TCGA-BRCA Gene Expression Project
# ==============================

# Install Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "DESeq2", "WGCNA", "survival", "survminer", 
                       "pheatmap", "ggplot2", "dplyr"))

# Load Libraries
library(TCGAbiolinks)
library(DESeq2)
library(WGCNA)
library(survival)
library(survminer)
library(pheatmap)
library(ggplot2)
library(dplyr)

# ------------------------------
# 1. Download and Prepare TCGA BRCA Data
# ------------------------------
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  sample.type = c("Primary Tumor", "Solid Tissue Normal"))
GDCdownload(query)
data <- GDCprepare(query)

# ------------------------------
# 2. Differential Expression Analysis
# ------------------------------
counts <- assay(data)
metadata <- colData(data)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ sample_type)
dds <- DESeq(dds)

res <- results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))
res <- lfcShrink(dds, coef="sample_type_Primary_Tumor_vs_Solid_Tissue_Normal", res=res)
sig_genes <- res[which(res$padj < 0.05), ]
write.csv(as.data.frame(sig_genes), "DEGs.csv")

# Volcano Plot
plot(res$log2FoldChange, -log10(res$pvalue), pch=20, 
     col=ifelse(res$padj < 0.05, "red", "black"),
     main="Volcano Plot")

# ------------------------------
# 3. Clustering and Heatmap
# ------------------------------
top_genes <- head(rownames(sig_genes[order(sig_genes$padj), ]), 50)
norm_counts <- counts(dds, normalized=TRUE)
pheatmap(norm_counts[top_genes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         main = "Top 50 DEGs Heatmap")

# ------------------------------
# 4. Co-expression Network Analysis (WGCNA)
# ------------------------------
options(stringsAsFactors = FALSE)
datExpr <- t(norm_counts[rowMeans(norm_counts) > 10, ])
datExpr <- datExpr[sample(1:nrow(datExpr), 100), ]  # Subsample for speed

powers <- c(1:10)
sft <- pickSoftThreshold(datExpr, powerVector = powers)
softPower <- sft$powerEstimate

net <- blockwiseModules(datExpr, power = softPower,
                        TOMType = "unsigned", minModuleSize = 30,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE)
plotDendroAndColors(net$dendrograms[[1]], net$colors, "Module Colors")

# ------------------------------
# 5. Survival Analysis
# ------------------------------
clinical <- GDCquery_clinic("TCGA-BRCA", type = "clinical")
gene <- "TP53"
expr <- norm_counts[gene, ]
sample_ids <- substr(names(expr), 1, 12)
clinical_expr <- data.frame(expr = expr, bcr_patient_barcode = sample_ids)

merged <- merge(clinical_expr, clinical, by = "bcr_patient_barcode")
merged <- merged[!is.na(merged$days_to_death), ]
merged$group <- ifelse(merged$expr > median(merged$expr), "High", "Low")

surv_obj <- Surv(time = as.numeric(merged$days_to_death), event = merged$vital_status == "Dead")
fit <- survfit(surv_obj ~ group, data = merged)

ggsurvplot(fit, data = merged, pval = TRUE, risk.table = TRUE, 
           title = paste("Survival based on", gene, "expression"))

# ------------------------------
# 6. Clinical and Expression Integration
# ------------------------------
expr_df <- data.frame(bcr_patient_barcode = sample_ids, expression = as.numeric(expr))
merged_data <- merge(expr_df, clinical, by = "bcr_patient_barcode")

# Expression vs Tumor Stage
ggplot(merged_data, aes(x = clinical_stage, y = expression)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = paste(gene, "Expression vs Tumor Stage"),
       x = "Tumor Stage", y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Expression vs Age
cor_test <- cor.test(as.numeric(merged_data$age_at_index), merged_data$expression, use = "complete.obs")
print(cor_test)

ggplot(merged_data, aes(x = age_at_index, y = expression)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", col = "red") +
  labs(title = paste(gene, "Expression vs Age"), x = "Age", y = "Expression Level")
