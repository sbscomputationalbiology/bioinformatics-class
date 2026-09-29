# 01_setup_renv.R
# Run line by line in the R console (not with source()), from the root of
# your project (setwd() there first).
#
# renv::init() restarts the R session, so the script is split in two parts:
#   PART 1: run once.
#   PART 2: run AFTER the restart, inside the renv project.

# ---- PART 1: create the renv project ----------------------------------------

# Install renv itself (outside the project library)
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# bioconductor = TRUE pins the project to the Bioconductor release that
# matches your R version. The R session restarts after this line.
renv::init(bioconductor = TRUE)


# ---- PART 2: install and verify packages (run after the restart) ------------

# Make sure the renv project is active in this session
if (is.null(renv::project())) {
  stop("No active renv project. setwd() to the project folder and restart R.")
}
.libPaths()   # first entry should be inside <project>/renv/library/

# Install the newest versions available for this Bioconductor release.
# The "bioc::" prefix makes renv pull them from Bioconductor.
#   Biostrings: reads/writes sequences (scripts 02 and 03)
#   pwalign:    pairwise alignment (script 03; split out of Biostrings in
#               recent Bioconductor releases)
renv::install(c("bioc::Biostrings", "bioc::pwalign"))

# Verify that they really installed. Stops with a clear message if not
# (a missing compiler is the usual cause; see notes at the bottom).
for (pkg in c("Biostrings", "pwalign")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " did not install. Read the error messages printed by ",
         "renv::install() above.")
  }
}

library(Biostrings)
packageVersion("Biostrings")
packageVersion("pwalign")
BiocManager::version()

# Record the exact versions in renv.lock
renv::snapshot()

# Notes: these packages and their dependencies (S4Vectors, IRanges,
# XVector, ...) contain compiled code. If no binary exists for your R
# version they are built from source, which needs a compiler:
#   Windows: Rtools matching your R version
#   macOS:   xcode-select --install
#   Linux:   build-essential (plus zlib1g-dev)
# Install the toolchain, then re-run PART 2.
# 

