#!/usr/bin/env julia

###############################################################################
#                                    README
# 
# This program read PDB structures and prepare toppology and coordinate files
# for CG MD simulations in Genesis.
#
# PDB format:
# 1. Atoms startswith "ATOM  "
# 2. Chains should end with "TER" and have different IDs
# 
# Unit in the script: kcal/mol, Å
# Unit for output:    kJ/mol, nm
###############################################################################

include("./lib/constants.jl")

using .GCGConstants

