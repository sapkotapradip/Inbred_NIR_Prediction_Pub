################################################################################
# AlphaSimR Breeding Program Simulation with Real Sorghum Genomic Data
# Sapkota, Pradip
#
# Description:
#   Simulates a two-heterotic-pool sorghum hybrid breeding program using:
#     - Real SNP genotype data (T92I012_clean.txt): 93 inbreds x ~38,234 markers
#     - Observed phenotypic BLUEs for grain yield, plant height, days to anthesis
#   Workflow:
#     1. Format real genotype/map data for AlphaSimR
#     2. Import inbred founders, split into A-pool (female) and R-pool (male)
#     3. Calibrate additive genetic variance from observed data
#     4. Simulate test-cross evaluation and selection over multiple cycles
#     5. Optionally layer genomic selection (GBLUP) using rrBLUP
#     6. Track and plot genetic gain, variance, and inbreeding per cycle
################################################################################

library(AlphaSimR)   # breeding simulation engine
library(dplyr)       # data wrangling
library(ggplot2)     # visualization

# ── 0. User-configurable parameters ──────────────────────────────────────────

set.seed(42)

DATA_DIR   <- file.path("..", "data")   # adjust if running from a different wd
N_CYCLES   <- 5                         # number of selection cycles to simulate
N_CROSSES  <- 20                        # random crosses within each pool per cycle
N_SELFINGS <- 5                         # SSD generations to advance new inbreds
N_TESTCROSS_TESTERS <- 5                # number of testers from opposite pool
SELECT_TOP  <- 0.20                     # fraction of hybrids selected as parents
USE_GS      <- TRUE                     # switch: genomic selection instead of PS

# ── 1. Load real genotype data ────────────────────────────────────────────────

cat("Loading SNP marker data...\n")
geno_raw <- read.table(
  file        = file.path(DATA_DIR, "T92I012_clean.txt"),
  header      = TRUE,
  check.names = FALSE
)

# Row names = line IDs; columns 2:end = markers coded 0/2 (homozygous inbreds)
line_ids  <- geno_raw[[1]]
geno_mat  <- as.matrix(geno_raw[, -1])   # lines x markers, values in {0, 2}
rownames(geno_mat) <- line_ids
marker_names <- colnames(geno_mat)

cat(sprintf("  %d inbred lines x %d markers loaded.\n",
            nrow(geno_mat), ncol(geno_mat)))

# ── 2. Build genetic map from marker names (format: "chrXX_<bp_position>") ──

# Sorghum has 10 chromosomes; assume ~1 cM per Mb recombination rate
BP_PER_MORGAN <- 1e8   # 100 Mb ≈ 1 Morgan (conservative estimate for sorghum)

cat("Building genetic map from marker positions...\n")
chr_labels <- sub("_[0-9]+$", "", marker_names)         # "chr01", "chr02", ...
bp_pos     <- as.numeric(sub("^chr[0-9]+_", "", marker_names))   # bp integer
chr_num    <- as.integer(sub("chr0?", "", chr_labels))  # 1 – 10

map_df <- data.frame(
  marker = marker_names,
  chr    = chr_num,
  bp     = bp_pos,
  stringsAsFactors = FALSE
)

# Convert bp to Morgans, reset to 0 at the start of each chromosome
map_df <- map_df %>%
  group_by(chr) %>%
  mutate(morgan = (bp - min(bp)) / BP_PER_MORGAN) %>%
  ungroup()

cat(sprintf("  Chromosomes: %s\n", paste(sort(unique(map_df$chr)), collapse = ", ")))
cat(sprintf("  Markers per chromosome: %s\n",
            paste(table(map_df$chr), collapse = ", ")))

# ── 3. Convert genotype coding: 0/2 → 0/1 haplotypes for AlphaSimR ───────────
# Inbreds are fully homozygous: allele 2 → 1, allele 0 → 0.
# AlphaSimR's newMapPop() expects haplotypes (one row per gamete; two rows per ind).

geno_01 <- geno_mat / 2L   # values now in {0, 1}

# Each inbred contributes two identical haplotypes (selfed to fixation)
n_lines   <- nrow(geno_01)
n_markers <- ncol(geno_01)

# Assemble haplotype matrix: 2*n_lines rows (odd = hap1, even = hap2)
hap_mat <- matrix(0L, nrow = 2 * n_lines, ncol = n_markers)
for (i in seq_len(n_lines)) {
  hap_mat[2 * i - 1, ] <- geno_01[i, ]
  hap_mat[2 * i,     ] <- geno_01[i, ]
}
rownames(hap_mat) <- rep(line_ids, each = 2)

# ── 4. Create AlphaSimR MapPop from real data ─────────────────────────────────

cat("Creating AlphaSimR MapPop from real genotypes...\n")

# Split haplotypes and map by chromosome
chr_list <- sort(unique(map_df$chr))

genMap  <- vector("list", length(chr_list))
hapList <- vector("list", length(chr_list))

for (k in seq_along(chr_list)) {
  ch  <- chr_list[[k]]
  idx <- which(map_df$chr == ch)
  genMap[[k]]  <- map_df$morgan[idx]
  hapList[[k]] <- hap_mat[, idx]
}

founderPop <- newMapPop(
  genMap   = genMap,
  haplotypes = hapList,
  inbred   = TRUE,    # founders are fully inbred (no residual heterozygosity)
  ploidy   = 2L
)

cat(sprintf("  MapPop created: %d individuals, %d total markers.\n",
            founderPop@nInd, sum(founderPop@nLoci)))

# ── 5. Define heterotic pools ─────────────────────────────────────────────────
# A-pool (female, cytoplasmic male-sterile maintainer lines): prefixes A*, AHF*, AKS*, AN*, ATx*
# R-pool (male, restorer lines): prefixes R*, RTx*

A_ids <- line_ids[grepl("^(A[0-9]+|AHF|AKS|AN[0-9]+|ATx)", line_ids)]
R_ids <- line_ids[grepl("^(R[0-9]+|RTx)",                  line_ids)]

cat(sprintf("  A-pool (female): %d lines | R-pool (male): %d lines\n",
            length(A_ids), length(R_ids)))

# Indices within the founder population (1-based)
A_idx <- which(line_ids %in% A_ids)
R_idx <- which(line_ids %in% R_ids)

poolA <- selectInd(founderPop, individuals = A_idx)
poolR <- selectInd(founderPop, individuals = R_idx)

# ── 6. Set up SimParam and trait architecture ─────────────────────────────────
# Calibrate additive variance from observed BLUEs for grain yield
# (average across 20CS and 20MA environments)

cat("Calibrating trait parameters from observed phenotypic data...\n")

pheno_20CS <- read.csv(file.path(DATA_DIR, "blues_20CS_pheno_parents.csv"))
pheno_20MA <- read.csv(file.path(DATA_DIR, "blues_20MA_pheno_parents.csv"))

calc_var <- function(df1, df2, col) {
  x <- c(df1[[col]], df2[[col]])
  x <- x[!is.na(x)]
  var(x)
}

varP_gy <- calc_var(pheno_20CS, pheno_20MA, "gy")   # phenotypic variance GY
varP_ph <- calc_var(pheno_20CS, pheno_20MA, "ph")   # phenotypic variance PH
varP_dy <- calc_var(pheno_20CS, pheno_20MA, "dy")   # phenotypic variance DA

# Assume broad-sense H2 ≈ 0.50 for GY, 0.65 for PH, 0.70 for DA (sorghum lit.)
H2_gy <- 0.50; H2_ph <- 0.65; H2_dy <- 0.70
varG_gy <- varP_gy * H2_gy
varG_ph <- varP_ph * H2_ph
varG_dy <- varP_dy * H2_dy
varE_gy <- varP_gy * (1 - H2_gy)
varE_ph <- varP_ph * (1 - H2_ph)
varE_dy <- varP_dy * (1 - H2_dy)

cat(sprintf("  GY: varG=%.1f, varE=%.1f\n", varG_gy, varE_gy))
cat(sprintf("  PH: varG=%.3f, varE=%.3f\n", varG_ph, varE_ph))
cat(sprintf("  DA: varG=%.3f, varE=%.3f\n", varG_dy, varE_dy))

SP <- SimParam$new(founderPop)

# Add three additive traits; nQtlPerChr = 50 QTL per chromosome
SP$addTraitA(
  nQtlPerChr = 50,
  mean       = c(mean(c(pheno_20CS$gy, pheno_20MA$gy), na.rm = TRUE),
                 mean(c(pheno_20CS$ph, pheno_20MA$ph), na.rm = TRUE),
                 mean(c(pheno_20CS$dy, pheno_20MA$dy), na.rm = TRUE)),
  var        = c(varG_gy, varG_ph, varG_dy),
  name       = c("GY", "PH", "DA")
)

# Residual (environmental) variance for phenotyping
SP$setVarE(varE = c(varE_gy, varE_ph, varE_dy))

# Enable genomic selection marker tracking
SP$addSnpChip(nSnpPerChr = 100)   # 1,000 SNP chip markers for GS (100 per chromosome × 10 chromosomes)

# Re-create pools under the SimParam
poolA <- newPop(founderPop[A_idx], simParam = SP)
poolR <- newPop(founderPop[R_idx], simParam = SP)

# ── 7. Helper functions ───────────────────────────────────────────────────────

#' Make testcross hybrids: each A-line x N_TESTCROSS_TESTERS R testers
make_testcrosses <- function(pool_A, pool_R, n_testers = N_TESTCROSS_TESTERS) {
  # Sample random testers from R-pool
  tester_idx <- sample(pool_R@nInd, min(n_testers, pool_R@nInd))
  testers    <- selectInd(pool_R, individuals = tester_idx, simParam = SP)
  crosses    <- hybridCross(pool_A, testers, simParam = SP)
  crosses
}

#' Select top hybrids based on GY (trait index 1) and return parent indices
select_parents_ps <- function(hybrids, top_frac = SELECT_TOP) {
  gy_vals <- pheno(hybrids, simParam = SP)[, "GY"]
  n_sel   <- max(1L, round(length(gy_vals) * top_frac))
  order(gy_vals, decreasing = TRUE)[seq_len(n_sel)]
}

#' Genomic selection: train GBLUP on hybrids, predict parents, return top indices
select_parents_gs <- function(pool, hybrids, top_frac = SELECT_TOP) {
  if (!requireNamespace("rrBLUP", quietly = TRUE)) {
    message("rrBLUP not installed; falling back to phenotypic selection.")
    return(select_parents_ps(hybrids, top_frac))
  }
  library(rrBLUP)
  # Marker matrix for hybrids (centered)
  M_hyb <- pullSnpGeno(hybrids, simParam = SP) - 1   # {-1, 0, 1}
  y_gy  <- pheno(hybrids, simParam = SP)[, "GY"]
  # Train model
  fit   <- mixed.solve(y = y_gy, Z = M_hyb)
  # Marker matrix for candidate parents
  M_par <- pullSnpGeno(pool, simParam = SP) - 1
  gebv  <- M_par %*% fit$u + as.numeric(fit$beta)
  n_sel <- max(1L, round(length(gebv) * top_frac))
  order(gebv, decreasing = TRUE)[seq_len(n_sel)]
}

# ── 8. Cycle-0 baseline: existing panel ───────────────────────────────────────

cat("\n--- Cycle 0 (founder panel) ---\n")

hyb0  <- make_testcrosses(poolA, poolR)
meanG <- c()
varG  <- c()

record_cycle <- function(pop, cycle_label) {
  gv <- gv(pop, simParam = SP)
  cat(sprintf("  %s | nInd=%d | mean(GY)=%.1f | var(GY)=%.1f\n",
              cycle_label, pop@nInd,
              mean(gv[, "GY"]),
              var(gv[, "GY"])
  ))
  list(cycle = cycle_label, nInd = pop@nInd,
       meanGY = mean(gv[, "GY"]), varGY  = var(gv[, "GY"]))
}

results <- list()
results[["C0"]] <- record_cycle(hyb0, "C0_hybrids")

# ── 9. Breeding cycles ────────────────────────────────────────────────────────

for (cyc in seq_len(N_CYCLES)) {

  cat(sprintf("\n--- Cycle %d ---\n", cyc))

  # (a) Make test-cross hybrids for evaluation
  hyb <- make_testcrosses(poolA, poolR)
  hyb <- setPheno(hyb, simParam = SP)   # add phenotypic noise

  # (b) Select superior parents
  if (USE_GS) {
    sel_A_idx <- select_parents_gs(poolA, hyb, SELECT_TOP)
    sel_R_idx <- select_parents_gs(poolR, hyb, SELECT_TOP)
  } else {
    # Phenotypic selection: evaluate parent pools directly and rank by GY
    poolA_pheno <- setPheno(poolA, simParam = SP)
    poolR_pheno <- setPheno(poolR, simParam = SP)
    n_sel_A <- max(1L, round(poolA@nInd * SELECT_TOP))
    n_sel_R <- max(1L, round(poolR@nInd * SELECT_TOP))
    sel_A_idx <- order(pheno(poolA_pheno, simParam = SP)[, "GY"],
                       decreasing = TRUE)[seq_len(n_sel_A)]
    sel_R_idx <- order(pheno(poolR_pheno, simParam = SP)[, "GY"],
                       decreasing = TRUE)[seq_len(n_sel_R)]
  }

  selA <- selectInd(poolA, individuals = sel_A_idx, simParam = SP)
  selR <- selectInd(poolR, individuals = sel_R_idx, simParam = SP)

  # (c) Generate new crosses within each pool
  newA_F1 <- randCross(selA, nCrosses = N_CROSSES, simParam = SP)
  newR_F1 <- randCross(selR, nCrosses = N_CROSSES, simParam = SP)

  # (d) Advance through SSD (successive self-fertilisation) to recover inbreds
  newA_inbred <- newA_F1
  newR_inbred <- newR_F1
  for (s in seq_len(N_SELFINGS)) {
    newA_inbred <- self(newA_inbred, simParam = SP)
    newR_inbred <- self(newR_inbred, simParam = SP)
  }

  # (e) Update pools: retain best existing lines + new inbreds (replace oldest)
  poolA <- c(selA, newA_inbred)
  poolR <- c(selR, newR_inbred)

  # (f) Record performance of hybrids this cycle
  results[[paste0("C", cyc)]] <- record_cycle(hyb, paste0("C", cyc, "_hybrids"))
}

# ── 10. Summarise and plot results ────────────────────────────────────────────

cat("\n=== Summary across cycles ===\n")
summary_df <- do.call(rbind, lapply(results, as.data.frame))
print(summary_df)

p_gain <- ggplot(summary_df, aes(x = cycle, y = meanGY, group = 1)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue") +
  labs(title  = "Genetic gain for grain yield across breeding cycles",
       x      = "Breeding cycle",
       y      = "Mean GY (simulated TBV)") +
  theme_bw(base_size = 13)

p_var <- ggplot(summary_df, aes(x = cycle, y = varGY, group = 1)) +
  geom_line(color = "tomato", linewidth = 1) +
  geom_point(size = 3, color = "tomato") +
  labs(title  = "Genetic variance for grain yield across breeding cycles",
       x      = "Breeding cycle",
       y      = "Var(GY) (simulated TBV)") +
  theme_bw(base_size = 13)

if (!dir.exists(file.path("..", "output"))) {
  dir.create(file.path("..", "output"), recursive = TRUE)
}
ggsave(file.path("..", "output", "AlphaSimR_genetic_gain.png"),  p_gain, width = 7, height = 4)
ggsave(file.path("..", "output", "AlphaSimR_genetic_variance.png"), p_var, width = 7, height = 4)

cat("\nPlots saved to output/.\n")

# ── 11. Optional: compare phenotypic vs genomic selection ─────────────────────

if (FALSE) {  # set to TRUE to run the comparison (doubles runtime)

  # Re-run cycles with phenotypic selection
  poolA_ps <- newPop(founderPop[A_idx], simParam = SP)
  poolR_ps <- newPop(founderPop[R_idx], simParam = SP)
  results_ps <- list()

  for (cyc in seq_len(N_CYCLES)) {
    hyb_ps       <- make_testcrosses(poolA_ps, poolR_ps)
    hyb_ps       <- setPheno(hyb_ps, simParam = SP)

    # Rank parents directly by their own phenotypic value (consistent with main loop)
    poolA_ps_ph  <- setPheno(poolA_ps, simParam = SP)
    poolR_ps_ph  <- setPheno(poolR_ps, simParam = SP)
    n_sel_A_ps   <- max(1L, round(poolA_ps@nInd * SELECT_TOP))
    n_sel_R_ps   <- max(1L, round(poolR_ps@nInd * SELECT_TOP))
    sel_A_ps     <- order(pheno(poolA_ps_ph, simParam = SP)[, "GY"],
                          decreasing = TRUE)[seq_len(n_sel_A_ps)]
    sel_R_ps     <- order(pheno(poolR_ps_ph, simParam = SP)[, "GY"],
                          decreasing = TRUE)[seq_len(n_sel_R_ps)]

    sA_ps        <- selectInd(poolA_ps, individuals = sel_A_ps, simParam = SP)
    sR_ps        <- selectInd(poolR_ps, individuals = sel_R_ps, simParam = SP)

    # Advance new inbreds through N_SELFINGS generations of selfing (same as main loop)
    newA_ps <- randCross(sA_ps, N_CROSSES, simParam = SP)
    newR_ps <- randCross(sR_ps, N_CROSSES, simParam = SP)
    for (s in seq_len(N_SELFINGS)) {
      newA_ps <- self(newA_ps, simParam = SP)
      newR_ps <- self(newR_ps, simParam = SP)
    }
    poolA_ps <- c(sA_ps, newA_ps)
    poolR_ps <- c(sR_ps, newR_ps)

    gv_ps      <- gv(hyb_ps, simParam = SP)
    results_ps[[paste0("C", cyc)]] <- data.frame(
      cycle  = paste0("C", cyc),
      meanGY = mean(gv_ps[, "GY"]),
      method = "PS"
    )
  }

  results_gs_df <- summary_df
  results_gs_df$method <- "GS"
  results_ps_df <- do.call(rbind, results_ps)
  compare_df    <- rbind(results_gs_df[, c("cycle", "meanGY", "method")],
                         results_ps_df)

  p_compare <- ggplot(compare_df, aes(x = cycle, y = meanGY,
                                       color = method, group = method)) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    labs(title = "GS vs. PS: genetic gain for grain yield",
         x = "Breeding cycle", y = "Mean GY (TBV)", color = "Method") +
    theme_bw(base_size = 13)

  ggsave(file.path("..", "output", "AlphaSimR_GS_vs_PS.png"),
         p_compare, width = 7, height = 4)
  cat("GS vs PS comparison plot saved.\n")
}
