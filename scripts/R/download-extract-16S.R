# 02_download_and_extract_16S.R
# Downloads the RefSeq assemblies for cases A-D from the NCBI FTP site into
#   <project folder>/raw-data/practice-8
# and runs the 16S extraction. Cases E and F are derived (simulated) from
# B/C and A, so they have no assembly to download.
#
# Portable: no absolute paths. The project folder is the renv project
# (e.g. .../bioinformatics-class), wherever it lives on each computer.
# Open that folder as the R project / working directory before running.

library(Biostrings)

options(timeout = 600)  # default 60 s can be too short on slow connections

# ---- Paths (relative to the project folder) ---------------------------------

raiz <- if (requireNamespace("renv", quietly = TRUE) && !is.null(renv::project())) {
  renv::project()
} else {
  getwd()
}

dir_datos <- file.path(raiz, "raw-data", "practice-8")
dir.create(dir_datos, recursive = TRUE, showWarnings = FALSE)
message("Data folder: ", normalizePath(dir_datos))

# ---- Genomes ----------------------------------------------------------------

genomas <- data.frame(
  caso      = c("A", "B", "C", "D"),
  organismo = c("Pseudomonas aeruginosa PAO1",
                "Escherichia coli K-12 MG1655",
                "Bacillus anthracis Ames",
                "Streptococcus pneumoniae TIGR4"),
  acc       = c("GCF_000006765.1", "GCF_000005845.2",
                "GCF_000007845.1", "GCF_000006885.1"),
  stringsAsFactors = FALSE
)

base_ftp <- "https://ftp.ncbi.nlm.nih.gov/genomes/all"

# GCF_000005845.2 -> https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845
ruta_acc <- function(acc) {
  num    <- sub("^GCF_([0-9]+)\\.[0-9]+$", "\\1", acc)
  partes <- substring(num, c(1, 4, 7), c(3, 6, 9))
  paste(c(base_ftp, "GCF", partes), collapse = "/")
}

# Look up the full assembly folder name (e.g. GCF_000005845.2_ASM584v2)
nombre_assembly <- function(acc) {
  html <- readLines(paste0(ruta_acc(acc), "/"), warn = FALSE)
  patron <- paste0(gsub(".", "\\.", acc, fixed = TRUE), "_[^/\"<]+")
  hits <- regmatches(html, regexpr(patron, html))
  if (length(hits) == 0) stop("Assembly folder not found for ", acc)
  unique(hits)[1]
}

for (i in seq_len(nrow(genomas))) {
  acc <- genomas$acc[i]
  asm <- nombre_assembly(acc)
  message("Caso ", genomas$caso[i], ": ", genomas$organismo[i], " (", asm, ")")
  
  # rna_from_genomic = rRNA/tRNA/ncRNA sequences (needed for the 16S code)
  # genomic          = full genome sequence
  for (sufijo in c("rna_from_genomic.fna.gz", "genomic.fna.gz")) {
    archivo <- paste0(asm, "_", sufijo)
    destino <- file.path(dir_datos, archivo)
    if (!file.exists(destino)) {
      download.file(paste(ruta_acc(acc), asm, archivo, sep = "/"),
                    destino, mode = "wb")
    }
  }
}

list.files(dir_datos)

# ---- 16S extraction ---------------------------------------------------------

extraer_16S <- function(archivo) {
  rna <- readDNAStringSet(archivo)
  rna[grepl("16S ribosomal RNA", names(rna))]
}

ecoli <- extraer_16S(file.path(dir_datos,
                               "GCF_000005845.2_ASM584v2_rna_from_genomic.fna.gz"))
length(ecoli)   # number of 16S copies in the genome
width(ecoli)    # length of each copy

# Save the first copy for submission to EzBioCloud (in the project folder,
# not in raw-data, which should only hold the downloaded files)
writeXStringSet(ecoli[1], file.path(raiz, "B_Ecoli_16S_copia1.fasta"))

# Copy number per genome, for all four cases
for (i in seq_len(nrow(genomas))) {
  f <- list.files(dir_datos,
                  pattern = paste0("^", genomas$acc[i], ".*rna_from_genomic"),
                  full.names = TRUE)
  x <- extraer_16S(f)
  message("Caso ", genomas$caso[i], ": ", length(x), " copies, widths ",
          paste(unique(width(x)), collapse = ", "))
}