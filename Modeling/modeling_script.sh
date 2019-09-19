#!/bin/bash

echo "Start: " `date`

# Execute Rosetta ab initio modeling with rosetta_options from rosetta_options file
~/software/ROSETTA/main/source/bin/AbinitioRelax.linuxgccrelease @rosetta_options

# Extract .pdb files from rosetta folding_silent files
~/software/ROSETTA/main/source/bin/extract_pdbs.mpi.linuxgccrelease -in:file:silent *out -in:file:fullatom -database ~/software/ROSETTA/main/database/

echo "Finish: " `date`
