# 04_length_and_region.R
# ==============================================================================
# Part D. Effect of length and variable region
#
# Goal: see how much the length and the region of the 16S gene matter for
# identification. Starting from ONE complete 16S copy (case B: E. coli, or
# case C: B. anthracis) we generate in silico the amplicons that two common
# primer pairs would give:
#   - 27F / 534R   -> V1-V3 region
#   - 515F / 806R  -> V4 region
# and save the complete gene plus both fragments as separate FASTA files, to
# be submitted one by one (by hand) to EzBioCloud. The results go in Table 2:
#   Sequence | Length (bp) | Top hit | % sim. | Compl. % |
#   No. of spp. >= 98.7 % | Correct identification?
# (Pay attention to completeness, number of species >= 98.7 % and any
# warning about length.)
#
# Requires: 01_setup_renv.R (Biostrings) and
#           02_download_and_extract_16S.R (files in raw-data/practice-8).
# Portable: every path is relative to the project folder.
# Output:   <project>/results/practice-8/EzBioCloud-3/*.fasta
# ==============================================================================

library(Biostrings)   # DNAString, matchPattern, subseq, FASTA writing

# ---- 1. Settings -------------------------------------------------------------

# Which genome to use: "B" (E. coli K-12 MG1655) or "C" (B. anthracis Ames)
caso <- "B"

casos <- list(
  B = list(acc = "GCF_000005845.2", nombre = "Ecoli"),
  C = list(acc = "GCF_000007845.1", nombre = "Banthracis")
)
genoma <- casos[[caso]]

# Primer pairs, written 5'->3' as usually published. Letters such as M, Y, N,
# V, W are IUPAC ambiguity codes (M = A/C, Y = C/T, N = any base, V = A/C/G,
# W = A/T).
cebadores <- list(
  V1V3 = list(fwd = "AGAGTTTGATCMTGGCTCAG", rev = "ATTACCGCGGCTGCTGG"),   # 27F / 534R
  V4   = list(fwd = "GTGYCAGCMGCCGCGGTAA",  rev = "GGACTACNVGGGTWTCTAAT") # 515F / 806R
)

# ---- 2. Paths (relative to the project folder) -------------------------------

# Project folder = active renv project; if there is none, the working directory
raiz <- if (requireNamespace("renv", quietly = TRUE) && !is.null(renv::project())) {
  renv::project()
} else {
  getwd()
}

# Input: genomes downloaded by script 02
dir_datos <- file.path(raiz, "raw-data", "practice-8")

# Output: FASTA files for the manual EzBioCloud upload (folder created if
# missing; this is the third upload folder, after EzBioCloud-2 of Part C).
# Change the folder name here if you want another one.
dir_salida <- file.path(raiz, "results", "practice-8", "EzBioCloud-3")
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

# Simulates a PCR on sequence s and returns the amplicon WITHOUT the primers.
#  - fwd: forward primer, matches the sequence as written.
#  - rev: reverse primer; it anneals to the opposite strand, so what we look
#    for in s is its reverse complement.
#  - fixed = FALSE lets the IUPAC ambiguity codes of the primers match.
#  - mm: mismatches allowed (0 = exact match apart from the ambiguity codes).
#    Raise it if a primer is "not found" in a divergent genome.
#  - The amplicon runs from the base after the forward primer to the base
#    before the reverse-primer site, i.e. primers excluded.
recortar <- function(s, fwd, rev, mm = 0) {
  f <- matchPattern(DNAString(fwd), s, max.mismatch = mm, fixed = FALSE)
  r <- matchPattern(reverseComplement(DNAString(rev)), s,
                    max.mismatch = mm, fixed = FALSE)
  if (length(f) == 0 || length(r) == 0) stop("Primer not found")
  
  # Only reverse-primer sites downstream of the first forward site make sense
  n_sitios <- c(fwd = length(f), rev = length(r))
  r <- r[start(r) > end(f)[1]]
  if (length(r) == 0) stop("Reverse primer not found downstream of the forward primer")
  if (any(n_sitios > 1)) {
    warning("A primer matches more than one site; the first valid one is used.")
  }
  
  message("  forward primer at ", start(f)[1], "-", end(f)[1],
          "; reverse primer site at ", start(r)[1], "-", end(r)[1])
  subseq(s, start = end(f)[1] + 1, end = start(r)[1] - 1)   # primers excluded
}

# ---- 4. Complete sequence ----------------------------------------------------

# All 16S copies of the chosen genome; we work with the first copy, so the
# "Full gene" row of Table 2 and both fragments come from the same sequence.
copias <- extraer_16S(archivo_rna(genoma$acc))
s <- copias[[1]]
message("Using 16S copy 1 of case ", caso, " (", names(copias)[1], ")")

# ---- 5. In silico fragments --------------------------------------------------

fragmentos <- lapply(names(cebadores), function(region) {
  message("Region ", region, ":")
  recortar(s, cebadores[[region]]$fwd, cebadores[[region]]$rev)
})
names(fragmentos) <- names(cebadores)

# Lengths in bp for the "Length" column of Table 2
longitudes <- c(complete = length(s), sapply(fragmentos, length))
longitudes

# ---- 6. Save FASTA files for EzBioCloud --------------------------------------

# One sequence per file (EzBioCloud is used one sequence at a time):
#   E_<organism>_completo.fasta, E_<organism>_V1V3.fasta, E_<organism>_V4.fasta
# "E" = derived sequences of case B or C (Case E of the practice).
completo <- DNAStringSet(list(completo = s))
guardar <- function(secuencia, etiqueta) {
  archivo <- file.path(dir_salida, paste0("E_", genoma$nombre, "_", etiqueta, ".fasta"))
  writeXStringSet(secuencia, archivo)
  message("Saved: ", archivo)
}

guardar(completo, "completo")
for (region in names(fragmentos)) {
  guardar(DNAStringSet(setNames(list(fragmentos[[region]]), region)), region)
}

# Both fragments together in a single multi-FASTA file (optional, for
# reference only: EzBioCloud needs them submitted separately)
frag <- DNAStringSet(fragmentos)
writeXStringSet(frag, file.path(dir_salida, paste0("E_", genoma$nombre, "_fragmentos.fasta")))

# ---- 7. Table 2 --------------------------------------------------------------

# After submitting each file to EzBioCloud, record for each row:
#   top hit, % similarity, completeness (%), number of species >= 98.7 %,
#   whether the identification is correct, and any warning about length.
# Full gene = E_<organism>_completo.fasta, V1-V3 = ..._V1V3.fasta,
# V4 = ..._V4.fasta