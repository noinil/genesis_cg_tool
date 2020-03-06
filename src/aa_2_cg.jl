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
include("./lib/selection.jl")
include("./lib/conformation.jl")
include("./lib/coarse_graining_subroutines.jl")
include("./lib/parsers.jl")
include("./lib/coarse_graining.jl")


function aa_2_cg(args)

    # -----------------
    # Parsing arguments
    # -----------------
    verbose            = get(args, "verbose", false)
    pdb_name           = get(args, "pdb", "")
    gen_pwmcos_itp     = get(args, "pwmcos", false)
    gen_pwmcos_ns_itp  = get(args, "pwmcos-ns", false)
    do_output_psf      = get(args, "psf", false)
    do_output_cgpdb    = get(args, "cgpdb", false)
    do_output_sequence = get(args, "show-sequence", false)
    ff_protein_name    = get(args, "force-field-protein", "AICG2+")
    ff_DNA_name        = get(args, "force-field-DNA", "3SPN.2C")
    ff_RNA_name        = get(args, "force-field-RNA", "HT")

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
        error("Wrong force field for DNA.")
    end
    if haskey(FF_RNA_DICT, ff_RNA_name)
        ff_rna = FF_RNA_DICT[ff_RNA_name]
    else
        error("Wrong force field for RNA.")
    end

    if gen_pwmcos_itp
        ff_pro_dna = FF_PWMcos
    elseif gen_pwmcos_ns_itp
        ff_pro_dna = FF_PWMcos_ns
    elseif get(args, "protein-DNA-Go", false)
        ff_pro_dna = FF_pro_DNA_Go
    else
        ff_pro_dna = FF_UNKNOWN
    end
    ff_pro_rna = FF_pro_RNA_Go
    ff_dna_rna = FF_UNKNOWN

    force_field = ForceFieldCG(ff_pro, ff_dna, ff_rna, ff_pro_dna, ff_pro_rna, ff_dna_rna)

    # -----------------------------------
    # FF modeling options from .toml file
    # -----------------------------------
    toml_filename = get(args, "config", "")
    if length(toml_filename) > 0
        new_toml_config = read_TOML(toml_filename)
        args["modeling-options"] = new_toml_config
    end


    ###########################################################################
    #                              Core Functions                             #
    ###########################################################################

    mol_name = get(args, "output-name", "")
    if length(mol_name) == 0
        mol_name = pdb_name[1:end-4] * "_cg"
    end

    # --------
    # Read PDB
    # --------
    if verbose
        println("============================================================")
        println("> Open PDB file:")
    end

    aa_molecule = read_PDB(pdb_name)

    aa_num_atom    = length(aa_molecule.atom_names)
    aa_num_residue = length(aa_molecule.residues)
    aa_num_chain   = length(aa_molecule.chains)

    if verbose
        println("          > Number of atoms    : $(aa_num_atom)")
        println("          > Number of residues : $(aa_num_residue)")
        println("          > Number of chains   : $(aa_num_chain)")
    end

    if do_output_sequence
        write_sequence(aa_molecule, mol_name)

        if verbose
            println("------------------------------------------------------------")
            println("[1;32m DONE! [0m ")
            println("============================================================")
        end

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

    if verbose
        println("============================================================")
        println("> Output CG .itp and .gro files.")
    end

    if do_output_top
        write_grotop(cg_top, mol_name, args)
    end
    if do_output_pwmcos
        write_grotop_pwmcos(cg_top, mol_name, args)
    end

    if do_output_gro
        write_grocrd(cg_top, cg_conf, mol_name, args)
    end

    if do_output_psf
        write_psf(cg_top, mol_name, args)
    end

    if do_output_cgpdb
        write_pdb(cg_top, cg_conf, mol_name, args)
    end

    if verbose
        println("------------------------------------------------------------")
        println("------------------------------------------------------------")
        println(">[1;32m FINISH! [0m ")
        println(" Please check the .itp and .gro files.")
        println("============================================================")
    end

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
        default = "HT"

        "--config", "-f"
        help = "Force field configuration details."
        arg_type = String
        default = ""

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

        "--use-safe-dihedral"
        help = "Safe dih potential: 0) do nothing (default); 1) remove dih w/ large angles; 2) sin(kθ) type; 3) sin^3(θ) type."
        arg_type = Int
        default = 0

        "--3spn-param"
        help = "Generate 3SPN.2C parameters from x3DNA generated PDB structure."
        arg_type = Int
        default = 0

        "--protein-DNA-Go"
        help = "Generate parameters for protein-DNA Go-like contact interactions."
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

        "--pwmcos-ns"
        help = "Generate parameters for protein-DNA sequence-NON-specific interactions."
        action = :store_true

        "--psf"
        help = "Prepare PSF file."
        action = :store_true

        "--cgpdb"
        help = "Prepare CG PDB file."
        action = :store_true

        "--cgconnect"
        help = "Prepare CG PDB file with CONECTed bonds."
        action = :store_true

        "--cgRNA-phosphate-Go"
        help = "Include phosphate in Go-type contact interactions."
        action = :store_true

        "--pfm", "-p"
        help = "Position frequency matrix file for protein-DNA sequence-specific interactions."
        arg_type = String
        default = ""

        "--test-local-only"
        help = "TEST: only generate local interaction parameters."
        action = :store_true

        "--patch"
        help = "Append (apply patch) to .itp file."
        arg_type = String
        default = ""

        "--show-sequence"
        help = "Show sequence of molecules in PDB."
        action = :store_true

        "--output-name"
        help = "Specify the system name for output."
        arg_type = String
        default = ""

        "--verbose", "-v"
        help = "Output more information."
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
