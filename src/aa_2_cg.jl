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

using Printf
using ArgParse
using Formatting

include("./lib/biomath.jl")
include("./lib/constants.jl")
include("./lib/molecule.jl")
include("./lib/topology.jl")
include("./lib/conformation.jl")
include("./lib/interactions.jl")
include("./lib/file_io.jl")
include("./lib/coarse_graining.jl")


function aa_2_cg(args)

    # -----------------
    # Parsing arguments
    # -----------------
    pdb_name           = args["pdb"]
    gen_pwmcos_itp     = args["pwmcos"]
    do_output_psf      = args["psf"]
    do_output_cgpdb    = args["cgpdb"]
    do_output_sequence = args["show-sequence"]
    ff_protein_name    = args["force-field-protein"]
    ff_DNA_name        = args["force-field-DNA"]
    ff_RNA_name        = args["force-field-RNA"]

    # ---------------
    # Set force field
    # ---------------
    if haskey(FF_PRO_DICT, ff_protein_name)
        ff_pro = FF_PRO_DICT[ff_protein_name]
    else
        error("Wrong force field for protein.")
    end
    if haskey(FF_DNA_DICT, ff_DNA_name)
        ff_dna = FF_DNA_DICT[ff_DNA_name]
    else
        error("Wrong force field for protein.")
    end
    if haskey(FF_RNA_DICT, ff_RNA_name)
        ff_rna = FF_RNA_DICT[ff_RNA_name]
    else
        error("Wrong force field for protein.")
    end

    if gen_pwmcos_itp
        ff_pro_dna = FF_PWMcos
    else
        ff_pro_dna = FF_UNKNOWN
    end
    force_field = ForceFieldCG(ff_pro, ff_dna, ff_rna, ff_pro_dna, FF_UNKNOWN, FF_UNKNOWN)


    ###########################################################################
    #                              Core Functions                             #
    ###########################################################################

    mol_name = pdb_name[1:end-4]

    # --------
    # Read PDB
    # --------
    println("============================================================")
    println("> Open PDB file:")

    aa_molecule  = read_PDB(pdb_name)

    aa_num_atom    = length(aa_molecule.atom_names)
    aa_num_residue = length(aa_molecule.residues)
    aa_num_chain   = length(aa_molecule.chains)

    println("          > Number of atoms    : $(aa_num_atom)")
    println("          > Number of residues : $(aa_num_residue)")
    println("          > Number of chains   : $(aa_num_chain)")

    if do_output_sequence
        write_sequence(aa_molecule, mol_name)
        return 0
    end

    # -----------------------------------
    # Make a CG topology from AA molecule
    # -----------------------------------
    cg_top, cg_conf = coarse_graining(aa_molecule, force_field)

    rg_all = radius_of_gyration(cg_conf)
    rc_all = radius_of_circumshpere(cg_conf)

    ###########################################################################
    #                            Output CG Topology                           #
    ###########################################################################

    if gen_pwmcos_itp
        do_output_top    = false
        do_output_itp    = false
        do_output_gro    = false
        do_output_pwmcos = true
    else
        do_output_top    = true
        do_output_itp    = true
        do_output_gro    = true
        do_output_pwmcos = false
    end

    println("============================================================")
    println("> Output CG .itp and .gro files.")


    println("------------------------------------------------------------")
    println("------------------------------------------------------------")
    println("[1;32m FINISH! [0m ")
    println(" Please check the .itp and .gro files.")
    println("============================================================")


end

