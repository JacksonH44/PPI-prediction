library(frustratometeR)
library(reticulate)

use_python(".venv/bin/python")
Sys.setenv(RETICULATE_PYTHON = ".venv/bin/python")
reticulate::py_config()

args = commandArgs(trailingOnly=TRUE)

PdbFile <- args[1]  # specify full path
ResultsDir <- args[2]  # specify full path
Chain <- args[3]  # PDB chain ("A" or "B")

Pdb_conf <- calculate_frustration(
    PdbFile = PdbFile, 
    Mode = "singleresidue", 
    ResultsDir = ResultsDir, 
    Chain = Chain, 
    Graphics = FALSE
)
