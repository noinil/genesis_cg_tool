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

    aa_molecule = read_aaPDB(pdb_name)

    aa_num_atom    = length(aa_molecule.atom_names)
    aa_num_residue = length(aa_molecule.residues)
    aa_num_chain   = length(aa_molecule.chains)

    println("          > Number of atoms    : $(aa_num_atom)")
    println("          > Number of residues : $(aa_num_residue)")
    println("          > Number of chains   : $(aa_num_chain)")

    if do_output_sequence
        write_sequence(aa_molecule, mol_name)

        println("------------------------------------------------------------")
        println("[1;32m DONE! [0m ")
        println("============================================================")

        return 0
    end

    # -----------------------------------
    # Make a CG topology from AA molecule
    # -----------------------------------
    cg_top, cg_conf = coarse_graining(aa_molecule, force_field, args)

    rg_all = radius_of_gyration(cg_conf)
    rc_all = radius_of_circumshpere(cg_conf)

    ###########################################################################
    #                            Output CG Topology                           #
    ###########################################################################

    if gen_pwmcos_itp
        do_output_top    = false
        do_output_gro    = false
        do_output_pwmcos = true
    else
        do_output_top    = true
        do_output_gro    = true
        do_output_pwmcos = false
    end


    if do_output_top
        write_cg_grotop(cg_top, force_field, mol_name, args)
    end
    if do_output_pwmcos
        write_cg_grotop_pwmcos(cg_top, force_field, mol_name, args)
    end

    if do_output_gro
        write_cg_grocrd(cg_top, cg_conf, mol_name, args)
    end

    if do_output_psf
        write_cg_psf(cg_top, mol_name, args)
    end

    if do_output_cgpdb
        write_cg_pdb(cg_top, cg_conf, mol_name, args)
    end


    println("============================================================")
    println("> Output CG .itp and .gro files.")


    println("------------------------------------------------------------")
    println("------------------------------------------------------------")
    printstyled(" FINISH! \n", bold=true, color=:green)
    println(" Please check the .itp and .gro files.")
    println("============================================================")


end

# =============================
# Parsing Commandline Arguments
# =============================
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin

        "pdb"
        help     = "PDB file name."
        required = true
        arg_type = String

        "--force-field-protein"
        help = "Force field for protein."
        arg_type = String
        default = "AICG2+"

        "--force-field-DNA"
        help = "Force field for DNA."
        arg_type = String
        default = "3SPN.2C"

        "--force-field-RNA"
        help = "Force field for RNA."
        arg_type = String
        default = "Go"

        "--CCGO-contact-scale"
        help = "Scaling native contact interaction coefficient."
        arg_type = Float64
        default = 1.0

        "--respac", "-c"
        help = "RESPAC protein charge distribution data."
        arg_type = String
        default = ""

        "--aicg-scale"
        help = "Scale AICG2+ local interactions: 0) average; 1) general (default)."
        arg_type = Int
        default = 1

        "--3spn-param"
        help = "Generate 3SPN.2C parameters from x3DNA generated PDB structure."
        action = :store_true

        "--pwmcos"
        help = "Generate parameters for protein-DNA sequence-specific interactions."
        action = :store_true

        "--pwmcos-scale"
        help = "Energy scaling factor for PWMcos."
        arg_type = Float64
        default = 1.0

        "--pwmcos-shift"
        help = "Energy shifting factor for PWMcos."
        arg_type = Float64
        default = 0.0

        "--psf"
        help = "Prepare PSF file."
        action = :store_true

        "--cgpdb"
        help = "Prepare CG PDB file."
        action = :store_true

        "--cgconect"
        help = "Prepare CG PDB file with CONECTed bonds."
        action = :store_true

        "--pfm", "-p"
        help = "Position frequency matrix file for protein-DNA sequence-specific interactions."
        arg_type = String
        default = ""

        "--patch"
        help = "Append (apply patch) to .itp file."
        arg_type = String
        default = ""

        "--show-sequence"
        help = "Show sequence of molecules in PDB."
        action = :store_true

        "--log"
        help = "Output information to log file."
        action = :store_true

        "--debug"
        help = "DEBUG."
        action = :store_true
    end

    return parse_args(s)
end

# ====
# Main
# ====

function main()

    args = parse_commandline()

    aa_2_cg(args)

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
