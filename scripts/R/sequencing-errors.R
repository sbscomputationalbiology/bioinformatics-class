# 05_sequencing_errors.R
# ==============================================================================
# Part E. Effect of sequencing errors
#
# Goal: simulate low-quality reads by introducing random substitutions into a
# 16S sequence of Case A (P. aeruginosa PAO1) at four error rates (0.5, 1, 2
# and 3 %), submit each mutant to EzBioCloud by hand, record the similarity to
# the top hit, and plot similarity against error rate with the two usual
# taxonomic thresholds:
#   - 98.7 % : below it, a sequence looks like a new SPECIES
#   - 94.5 % : below it, a sequence looks like a new GENUS
#
# Requires: 01_setup_renv.R (Biostrings) and
#           02_download_and_extract_16S.R (files in raw-data/practice-8).
# Portable: every path is relative to the project folder.
# Output:   <project>/results/practice-8/EzBioCloud-4/*.fasta   (to upload)
#           <project>/results/practice-8/EzBioCloud-4/F_similitud_vs_error.png (plot)
# ==============================================================================

library(Biostrings)   # DNAString, DNAStringSet, FASTA reading/writing

# ---- 1. Settings -------------------------------------------------------------

semilla <- 2026                    # random seed: use your team number
tasas   <- c(0.005, 0.01, 0.02, 0.03)   # error rates: 0.5, 1, 2 and 3 %

# Similarity (%) to the top hit reported by EzBioCloud for each mutant, in the
# same order as 'tasas'. Leave NA until you have the results, then fill them in
# and run section 6 again.
sim_ezbiocloud <- c(NA, NA, NA, NA)

umbral_especie <- 98.7   # % similarity: species threshold
umbral_genero  <- 94.5   # % similarity: genus threshold

# ---- 2. Paths (relative to the project folder) -------------------------------

# Project folder = active renv project; if there is none, the working directory
raiz <- if (requireNamespace("renv", quietly = TRUE) && !is.null(renv::project())) {
  renv::project()
} else {
  getwd()
}

# Input: genomes downloaded by script 02
dir_datos <- file.path(raiz, "raw-data", "practice-8")

# Output: FASTA files for the manual EzBioCloud upload (fourth upload folder,
# after EzBioCloud-2 and EzBioCloud-3); the plot is saved in the same folder
dir_res    <- file.path(raiz, "results", "practice-8")
dir_salida <- file.path(dir_res, "EzBioCloud-4")
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)
message("Output folder: ", normalizePath(dir_salida))

# ---- 3. Helper functions -----------------------------------------------------

# Reads a RefSeq "rna_from_genomic" FASTA and keeps only the 16S rRNA records
extraer_16S <- function(archivo) {
  rna <- readDNAStringSet(archivo)
  rna[grepl("16S ribosomal RNA", names(rna))]
}

# Finds the rna_from_genomic file of an assembly accession in raw-data/practice-8
archivo_rna <- function(acc) {
  f <- list.files(dir_datos,
                  pattern = paste0("^", acc, ".*rna_from_genomic\\.fna\\.gz$"),
                  full.names = TRUE)
  if (length(f) != 1) {
    stop("Expected one rna_from_genomic file for ", acc, " in ", dir_datos,
         ". Run 02_download_and_extract_16S.R first.")
  }
  f
}

# Introduces random substitutions into sequence s at a given error rate.
#  - x: the sequence as a vector of single letters.
#  - pos: positions to mutate, drawn WITHOUT replacement, so exactly
#    round(length * tasa) different positions change.
#  - Each chosen base is replaced by one of the other three bases (never by
#    itself), so every draw is a real substitution.
# Only substitutions are simulated (no insertions or deletions).
mutar <- function(s, tasa) {
  x <- strsplit(as.character(s), "")[[1]]
  pos <- sample(length(x), round(length(x) * tasa))
  for (p in pos) x[p] <- sample(setdiff(c("A", "C", "G", "T"), x[p]), 1)
  DNAString(paste(x, collapse = ""))
}

# Number of positions at which two sequences of equal length differ
n_diferencias <- function(a, b) {
  sum(strsplit(as.character(a), "")[[1]] != strsplit(as.character(b), "")[[1]])
}

# ---- 4. Original sequence and mutants ----------------------------------------

# Fixing the seed makes the random substitutions reproducible
set.seed(semilla)

# First 16S copy of P. aeruginosa PAO1 (Case A)
pae <- extraer_16S(archivo_rna("GCF_000006765.1"))[[1]]
length(pae)

# One mutant per error rate, named e.g. "Pae_error_0.5pct"
mut <- DNAStringSet(lapply(tasas, function(t) mutar(pae, t)))
names(mut) <- paste0("Pae_error_", tasas * 100, "pct")

# Sanity check: how many bases changed, and the exact identity of each mutant
# with the original sequence. The identity with the ORIGINAL is
# 100 * (1 - substitutions / length): the similarity EzBioCloud reports will be
# close to this value, because the top hit is almost identical to the original.
sustituciones <- sapply(mut, function(m) n_diferencias(pae, m))
tabla <- data.frame(
  error_pct        = tasas * 100,
  sustituciones    = sustituciones,
  identidad_pct    = round(100 * (1 - sustituciones / length(pae)), 2),
  sim_ezbiocloud   = sim_ezbiocloud,
  row.names        = names(mut)
)
tabla

# ---- 5. Save FASTA files for EzBioCloud --------------------------------------

# One sequence per file (submitted one at a time), e.g.
# F_Pae_error_0_5pct.fasta. Dots in the name are replaced by "_".
for (k in seq_along(mut)) {
  archivo <- file.path(dir_salida,
                       paste0("F_", gsub("\\.", "_", names(mut)[k]), ".fasta"))
  writeXStringSet(mut[k], archivo)
  message("Saved: ", archivo)
}

# All mutants together (optional, for reference)
writeXStringSet(mut, file.path(dir_salida, "F_mutantes.fasta"))

# ---- 6. Plot: similarity vs error rate ---------------------------------------

# Grey line = expected similarity if the top hit were identical to the
# original: similarity = 100 - error rate. Points = what EzBioCloud reports.
# Dashed / dotted horizontal lines = species / genus thresholds.
dibujar <- function() {
  plot(NA, xlim = c(0, 6), ylim = c(93, 100.5),
       xlab = "Simulated error rate (%)",
       ylab = "Similarity to top hit (%)",
       main = "Effect of sequencing errors on identification (Pae PAO1)")
  abline(a = 100, b = -1, col = "grey60")                   # expected similarity
  abline(h = umbral_especie, lty = 2, col = "red",  lwd = 2)  # species
  abline(h = umbral_genero,  lty = 3, col = "blue", lwd = 2)  # genus
  text(6, umbral_especie + 0.3, "98.7 % (species)", adj = 1, col = "red")
  text(6, umbral_genero  + 0.3, "94.5 % (genus)",   adj = 1, col = "blue")
  points(tabla$error_pct, tabla$identidad_pct, pch = 1, cex = 1.3, col = "grey40")
  points(tabla$error_pct, tabla$sim_ezbiocloud, pch = 19, cex = 1.3)
  legend("bottomleft", bty = "n",
         legend = c("Expected (100 - error)", "Exact identity with original",
                    "EzBioCloud top hit"),
         lty = c(1, NA, NA), pch = c(NA, 1, 19), col = c("grey60", "grey40", "black"))
}

dibujar()                                   # show it in the R plot window

png(file.path(dir_salida, "F_similitud_vs_error.png"), width = 1800, height = 1200, res = 200)
dibujar()                                   # save the same plot as a PNG
dev.off()

# ---- 7. Error rate at which the thresholds are crossed -----------------------

# With similarity = 100 - error rate, a known isolate falls below a threshold
# when the error rate exceeds 100 - threshold.
100 - umbral_especie   # error rate (%) above which it looks like a new species
100 - umbral_genero    # error rate (%) above which it looks like a new genus

