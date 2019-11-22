#!/usr/bin/env julia

using Random
using Printf
using ArgParse

include("../../../src/lib/biomath.jl")
include("../../../src/lib/molecule.jl")
include("../../../src/lib/topology.jl")
include("../../../src/lib/constants.jl")
include("../../../src/lib/selection.jl")
include("../../../src/lib/interactions.jl")
include("../../../src/lib/conformation.jl")
include("../../../src/lib/coarse_graining.jl")
include("../../../src/lib/parsers.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin

        "--sequence", "-s"
        help = "Protein sequence file."
        arg_type = String
        default = ""

        "--length"
        help = "Number of amino acids in the protein with random-sequence."
        arg_type = Int
        default = 100

        "--straightness"
        help = "Angle threshold to create a (non)-straight chain."
        arg_type = Float64
        default = 45.0

        "--strategy"
        help = "Strategy to construct the conformation of protein."
        arg_type = String
        default = "straight"

        "--force-field-protein"
        help = "Force field for protein."
        arg_type = String
        default = "AICG2+"

    end

    return parse_args(s)
end


function make_cg_protein_structure(args)

    seq_name = get(args, "sequence", "")

    # protein model
    ff_protein_name = get(args, "force-field-protein", "")
    if haskey(FF_PRO_DICT, ff_protein_name)
        ff_pro = FF_PRO_DICT[ff_protein_name]
    else
        error("Wrong force field for protein.")
    end

    # non-straightness (because ideal straight chain could have problem...)
    threshold_angle = get( args, "straightness", 45.0)

    println("============================================================")

    AA_FULLNAME_DICT = Dict(
        'A' => "ALA",
        'R' => "ARG",
        'N' => "ASN",
        'D' => "ASP",
        'C' => "CYS",
        'Q' => "GLN",
        'E' => "GLU",
        'G' => "GLY",
        'H' => "HIS",
        'I' => "ILE",
        'L' => "LEU",
        'K' => "LYS",
        'M' => "MET",
        'F' => "PHE",
        'P' => "PRO",
        'S' => "SER",
        'T' => "THR",
        'W' => "TRP",
        'Y' => "TYR",
        'V' => "VAL"
    )

    # ========================================================
    # Protein sequence (read from file or generate random one)
    # ========================================================
    if length(seq_name) > 0
        println("> Open sequence file:", seq_name)

        mol_name = split(basename( seq_name ), '.')[1]

        # -----------------------------------
        # read in protein sequence from fasta
        # -----------------------------------
        protein_seqence = ""
        num_chain = 0
        for line in eachline(seq_name)
            if length(line) == 0
                continue
            end
            if line[1] == '>'
                num_chain += 1
                continue
            end
            if num_chain > 1
                error("Only support single-chain protein!")
            end
            seq = strip(line)
            if length(seq) == 0
                continue
            end
            for b in seq
                if ! haskey(AA_FULLNAME_DICT, b)
                    error("Wrong protein sequence!")
                end
            end
            protein_seqence *= seq
        end
        protein_length = length(protein_seqence)
    else
        println("> Generating random protein sequence:")

        mol_name = "random_protein"

        # --------------------------------
        # generate random protein sequence
        # --------------------------------
        protein_length = args["length"]
        protein_seqence = randstring("ARNDCQEGHILKMFPSTWYVX", protein_length)
    end
    println("> Protein sequence:   ( Length: $protein_length ) ")
    println("> ", protein_seqence)

    # =======================================
    # Generating AAMolecule based on sequence
    # =======================================
    atom_names = Vector{String}(undef, protein_length)
    atom_coors = zeros(3, protein_length)
    residues   = Vector{AAResidue}(undef, protein_length)
    chains     = Vector{AAChain}(undef, 1)

    if args["strategy"] == "straight"
        for i in 1 : protein_length
            aa_short_name = protein_seqence[i]
            aa_residue_name = AA_FULLNAME_DICT[aa_short_name]
            
            atom_names[i] = "CA"
            if i > 1
                theta = rand() * threshold_angle
                phi   = ( rand() - 0.5 ) * threshold_angle
                atom_coors[:, i] = atom_coors[:, i - 1] + [cosd(theta), sind(theta) * cosd(phi), sind(theta) * sind(phi)] * 3.8
            end
            residues[i] = AAResidue(aa_residue_name, [i])
        end
        new_chain = AAChain('A', rpad(mol_name, 4)[1:4], MOL_PROTEIN, [i for i = 1 : protein_length])
        chains[1] = new_chain
    elseif args["strategy"] == "random-walk"
        println("Self-avoiding random walk not support yet.")
    end

    new_mol = AAMolecule(atom_names, atom_coors, residues, chains)

    # ===============================
    # coarse graining from AAMolecule
    # ===============================
    force_field = ForceFieldCG(ff_pro, 1, 1, 0, 0, 0)
    cg_top, cg_conf = coarse_graining(new_mol, force_field, args)

    args["modeling-options"] = Dict("IDR" => Dict("HPS_region" => "1 to $protein_length"))
    args["cgconect"] = true
    write_cg_grotop(cg_top, force_field, mol_name, args)
    write_cg_grocrd(cg_top, cg_conf, mol_name, args)
    write_cg_pdb(cg_top, cg_conf, mol_name, args)
    write_cg_psf(cg_top, mol_name, args)

    return 0
end


function main()
    args = parse_commandline()
    make_cg_protein_structure(args)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
