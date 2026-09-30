# 03_intragenomic_heterogeneity.R
# ==============================================================================
# Part C. Intragenomic heterogeneity
#
# Goal: bacterial genomes often carry several copies of the 16S rRNA gene
# (rrn operons), and those copies are not always identical. Here we
#   1. extract all 16S copies of a genome,
#   2. compare every copy with every other copy (pairwise global alignment),
#   3. find the two copies that differ the most, and
#   4. save those two copies as FASTA files, ready to be uploaded by hand to
#      EzBioCloud, to see whether they are identified as the same species.
#
# Cases: B (E. coli K-12 MG1655) and, if time permits, C (B. anthracis Ames).
#
# Requires: 01_setup_renv.R (Biostrings + pwalign) and
#           02_download_and_extract_16S.R (files in raw-data/practice-8).
# Portable: every path is relative to the project folder.
# Output:   <project>/results/practice-8/EzBioCloud-2/*.fasta
# ==============================================================================

library(Biostrings)   # reading/writing sequences (DNAStringSet, FASTA files)

# pwalign holds the alignment functions (pairwiseAlignment, pid). Stop early
# with a clear message if it is not installed in the renv project.
if (!requireNamespace("pwalign", quietly = TRUE)) {
  stop("Package 'pwalign' is missing. Run: renv::install(\"bioc::pwalign\"); renv::snapshot()")
}

# ---- 1. Paths (relative to the project folder) -------------------------------

# Project folder = active renv project; if there is none, the working directory.
raiz <- if (requireNamespace("renv", quietly = TRUE) && !is.null(renv::project())) {
  renv::project()
} else {
  getwd()
}

# Input: genomes downloaded by script 02
dir_datos <- file.path(raiz, "raw-data", "practice-8")

# Output: FASTA files for the manual EzBioCloud upload. recursive = TRUE
# creates results/, results/practice-8/ and results/practice-8/EzBioCloud-2/
# in one go if they do not exist yet.
dir_salida <- file.path(raiz, "results", "practice-8", "EzBioCloud-2")
dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)
message("Output folder: ", normalizePath(dir_salida))

# ---- 2. Helper functions -----------------------------------------------------

# Reads a RefSeq "rna_from_genomic" FASTA and keeps only the 16S rRNA records.
# The record headers contain "product=16S ribosomal RNA" for every copy.
extraer_16S <- function(archivo) {
  rna <- readDNAStringSet(archivo)
  rna[grepl("16S ribosomal RNA", names(rna))]
}

# Finds the rna_from_genomic file of an assembly accession (e.g.
# "GCF_000005845.2") inside raw-data/practice-8, whatever the assembly name.
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

# Builds an n x n matrix with the % identity between every pair of copies.
#  - type = "global": align the full length of both sequences (end to end).
#  - PID1: percent identity = identical positions / (aligned positions +
#    internal gap positions) x 100, so gaps also count as differences.
#  - The diagonal is 100 because each copy is compared with itself.
# Rows/columns are named c1..cn, following the order of the sequences.
matriz_pid <- function(secs) {
  n <- length(secs)
  etiquetas <- paste0("c", seq_len(n))
  M <- matrix(NA_real_, n, n, dimnames = list(etiquetas, etiquetas))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    aln <- pwalign::pairwiseAlignment(secs[[i]], secs[[j]], type = "global")
    M[i, j] <- pwalign::pid(aln, type = "PID1")
  }
  M
}

# Finds the pair of copies with the LOWEST identity in the matrix M.
#  - The matrix is symmetric, so we blank the diagonal and the lower triangle
#    (set to NA) to look at each pair only once and to ignore the trivial
#    self-comparisons.
#  - If several pairs tie for the minimum, the first one (lowest indices) is
#    used; n_empates reports how many pairs share that value.
copias_mas_divergentes <- function(M) {
  if (nrow(M) < 2) stop("At least two copies are needed to compare them.")
  M_sup <- M
  M_sup[lower.tri(M_sup, diag = TRUE)] <- NA
  minimo <- min(M_sup, na.rm = TRUE)
  pares <- which(M_sup == minimo, arr.ind = TRUE)
  pares <- pares[order(pares[, "row"], pares[, "col"]), , drop = FALSE]
  list(i = unname(pares[1, "row"]),
       j = unname(pares[1, "col"]),
       pid = minimo,
       n_empates = nrow(pares))
}

# Reports the most divergent pair and writes each of the two copies to its own
# FASTA file (one sequence per file, as needed for a manual EzBioCloud
# submission). File names: <prefijo>_16S_c<k>.fasta, where k is the copy
# number used in the matrix (c1, c2, ...). The original FASTA header is kept,
# so the locus tag / gene name identifies the copy in the genome.
guardar_pareja_divergente <- function(secs, M, prefijo) {
  p <- copias_mas_divergentes(M)
  
  message("Most divergent pair: c", p$i, " vs c", p$j,
          " (", round(p$pid, 2), " % identity)")
  if (p$n_empates > 1) {
    message("  Note: ", p$n_empates, " pairs share this value; the first one is used.")
  }
  if (p$pid == 100) {
    warning("All copies are identical: the 'most divergent' pair is arbitrary.")
  }
  message("  c", p$i, ": ", names(secs)[p$i], " (", width(secs)[p$i], " bp)")
  message("  c", p$j, ": ", names(secs)[p$j], " (", width(secs)[p$j], " bp)")
  
  for (k in c(p$i, p$j)) {
    archivo <- file.path(dir_salida, paste0(prefijo, "_16S_c", k, ".fasta"))
    writeXStringSet(secs[k], archivo)
    message("  Saved: ", archivo)
  }
  invisible(p)
}

# ---- 3. Case B: E. coli K-12 MG1655 -----------------------------------------

# Load the 16S copies of E. coli (7 are expected)
ecoli <- extraer_16S(archivo_rna("GCF_000005845.2"))
length(ecoli)

# All-vs-all identity matrix (values in %, rounded for display)
M_ecoli <- matriz_pid(ecoli)
round(M_ecoli, 2)

# Lowest identity between any two copies of the same genome, i.e. the greatest
# divergence within E. coli
min(M_ecoli)

# Identify the two most divergent copies and save them for EzBioCloud
guardar_pareja_divergente(ecoli, M_ecoli, prefijo = "B_Ecoli")

# ---- 4. Case C (if time permits): B. anthracis Ames -------------------------

# Same analysis for B. anthracis; delete this section if you are short on time
anthracis <- extraer_16S(archivo_rna("GCF_000007845.1"))
length(anthracis)   # number of 16S copies in this genome

M_anthracis <- matriz_pid(anthracis)
round(M_anthracis, 2)
min(M_anthracis)

guardar_pareja_divergente(anthracis, M_anthracis, prefijo = "C_Banthracis")