# Find differentially expressed genes with DESeq2 in Positron ----
# Author: Sofia Fertuzinhos, PhD
# Date: 2026

# ======================================== #
#                 PART I ----
# ======================================== #

# 1. Setup working directory ----
# Open the workshop folder (differential-expression) in Positron:
#   File -> Open Folder... -> select the differential-expression folder
# Positron starts R with that folder as the working directory. Check with:
#   getwd()
# here() builds every path in this script relative to that folder.
# Installation (run once): install.packages("here")

library(here)

# Results tables and figures are written to results/ and figures/, respectively 
# (created if it doesn't exist)
for (d in c("results", "figures")) {
  dir.create(here(d), showWarnings = FALSE)
}

# 2. Input data ----

# 2.1 Read raw counts ----
rawc <- read.csv(here("data", "clean", "shank3_rawcounts_clean.csv"),
                 header = TRUE,
                 row.names = 1
                 )

# 2.2 Read sample metadata ----
info <- read.csv(here("data", "clean", "shank3_metadata_clean.csv"),
                 header = TRUE,
                 stringsAsFactors = TRUE
                 )




# 3. Define reference levels for comparisons ----

# Establish the reference groups (baselines) used in differential expression tests.

# Comparisons are driven by two metadata columns:
#   - genotype
#   - Condition

# Use relevel() to set the reference level for each factor.

info$genotype  <- relevel(info$genotype, " ")
info$Condition <- relevel(info$Condition, " ")

# 4. Load DESeq2 (Bioconductor) ----
# DESeq2 is a Bioconductor package.
# Installation (run once):
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")

library(DESeq2)

# 5. Create the DESeq2 dataset object (DESeqDataSet) ----
# Build a DESeqDataSet from:
#   - rawc: counts matrix (genes x samples)
#   - info: sample metadata (samples x variables)
# The design formula specifies the variables used in the model.

Shank3 <- DESeqDataSetFromMatrix(
  countData =   ,
  colData   =   ,
  design    = ~ 
)


# ======================================== #
#                 PART II ----
# ======================================== #


## Ted Ellsworth presentation on  Design Matrix ----

Shank3 <- DESeqDataSetFromMatrix(
  countData = rawc,
  colData = info,
  design = ~ genotype + Condition + genotype:Condition
)


# ======================================== #
#                 PART III ----
# ======================================== #

# 6. Run DESeq2 differential expression analysis ----
# Run the DESeq2 pipeline (normalization, dispersion estimation, model fitting).
DEX_Shank3 <- DESeq(Shank3)

# 7. Quality control: PCA on transformed counts ----
# Use rlog() (or vst()) for visualization and sample-level QC.
log_DEX_Shank3 <- rlog(DEX_Shank3)

# Choose the metadata column(s) to color/group samples by in the PCA.
# Example: "genotype" or "Condition" or "Geno_Cond" (if you created that column).
plotPCA(log_DEX_Shank3, intgroup = "Condition")

# 8. Update the model design to include an interaction (if needed) ----
# If QC suggests an interaction between genotype and condition, update the design.
design(DEX_Shank3) <- ~ Condition + genotype + Condition:genotype

# 9. Inspect the model matrix (direction of comparisons) ----
# In the model matrix, 0 indicates the reference level and 1 the test level.
mod_mat <- stats::model.matrix(design(DEX_Shank3), colData(DEX_Shank3))
View(mod_mat)

# 10. List available DESeq2 result coefficients ----
# These names define what you can request via results(..., name=...) or contrast=...
resultsNames(DEX_Shank3)

#[1] "Intercept"                              
#[2] "genotype_Shank3_vs_WT"------------------- main genotype effect at the reference Condition level  
#[3] "Condition__sleep_deprived_vs__normal"---- main condition effect at reference genotype
#[4] "genotypeShank3.Condition_sleep_deprived"- interaction term

# 11. Extract DESeq2 results for specific comparisons ----

# 11.1 Default results() output (depends on the last term in the design) ----
res1 <- results(DEX_Shank3)
summary(res1)

# 11.2 Interaction term: genotype x condition ----
res2 <- results(
  DEX_Shank3,
  contrast = list("genotypeShank3.Condition_sleep_deprived")
)
summary(res2)

# 11.3 Shank3 vs WT (main genotype effect at the reference Condition level) ----
res3 <- results(
  DEX_Shank3,
  contrast = list("genotype_Shank3_vs_WT")
)
summary(res3)

# 11.4 Shank3 vs WT in sleep-deprived condition (main effect + interaction) ----
res4 <- results(
  DEX_Shank3,
  contrast = list(
    "genotype_Shank3_vs_WT",
    "genotypeShank3.Condition_sleep_deprived"
  )
)
summary(res4)

# 11.5 Sleep-deprived vs normal in WT (main Condition effect at reference genotype) ----
res5 <- results(
  DEX_Shank3,
  contrast = list("Condition__sleep_deprived_vs__normal")
)
summary(res5)

# 11.6 Sleep-deprived vs normal in Shank3 (main effect + interaction) ----
res6 <- results(
  DEX_Shank3,
  contrast = list(
    "Condition__sleep_deprived_vs__normal",
    "genotypeShank3.Condition_sleep_deprived"
  )
)
summary(res6)

# 12. Export results tables (.csv) ----
write.csv(res2, file = here("results", "Shank3_genoCond_interaction.csv"))
write.csv(res3, file = here("results", "Shank3_geno_Shank3vsWT_normal.csv"))
write.csv(res4, file = here("results", "Shank3_geno_Shank3vsWT_sleepDep.csv"))
write.csv(res5, file = here("results", "Shank3_Cond_SleepDepvsNormal_WT.csv"))
write.csv(res6, file = here("results", "Shank3_Cond_SleepDepvsNormal_Shank3.csv"))

# 13. Additional DESeq2 plots ----

log_DEX_Shank3 <- rlog(DEX_Shank3)
plotPCA(log_DEX_Shank3, intgroup = "Geno_Cond")


# 13.1 MA plot: log2 fold change vs mean normalized counts ----
Shank3_MA <- plotMA(
  DEX_Shank3,
  contrast = list("Condition__sleep_deprived_vs__normal")
)

# 13.2 Plot normalized counts for a single gene across groups ----
plotCounts(
  DEX_Shank3,
  gene     = "ENSMUSG00000002831|Plin4",
  intgroup = "Geno_Cond"
)

# 14. Export normalized counts ----
# Normalized counts = counts divided by sample-specific size factors.
normCounts <- counts(DEX_Shank3, normalized = TRUE)
write.csv(normCounts, file = here("results", "Shank3_normCounts_Deseq2.csv"))

# 15. Record session information for reproducibility ----
sessionInfo()
