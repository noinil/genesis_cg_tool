#!/usr/bin/env julia

using Random
using Printf
using ArgParse

include("../../../src/lib/gcj.jl")

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

        "--IDR-model"
        help = "IDR-model: 1) HPS; 2) KH; 3) AICG2+ LOCAL"
        arg_type = Int
        default = 1

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

    idr_model = get(args, "IDR-model", 1)

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

    chain_id_lib = "_ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

    # ========================================================
    # Protein sequence (read from file or generate random one)
    # ========================================================
    if length(seq_name) > 0
        println("> Open sequence file:", seq_name)

        mol_name = split(basename( seq_name ), '.')[1] * "_cg"

        # -----------------------------------
        # read in protein sequence from fasta
        # -----------------------------------
        num_chain, seq_list = read_fasta(seq_name)
        protein_seqence = seq_list[1]
        protein_length  = length(protein_seqence)

    else
        println("> Generating random protein sequence:")

        mol_name = "random_protein_cg"

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
    chain0     = Vector{AAChain}(undef, 1)

    if args["strategy"] == "straight"
        theta = rand() * threshold_angle
        for i in 1 : protein_length
            aa_short_name = protein_seqence[i]
            aa_residue_name = AA_FULLNAME_DICT[aa_short_name]

            atom_names[i] = "CA"
            if i == 1
                atom_coors[:, i] = [0, 0, 0]
            else
                if mod(i, 2) == 0
                    theta = rand() * threshold_angle
                    phi   = rand() * 360
                    atom_coors[:, i] = atom_coors[:, i - 1] + [sind(theta) * cosd(phi), sind(theta) * sind(phi), cosd(theta)] * 3.8
                else
                    atom_coors[:, i] = atom_coors[:, i - 2] + [0, 0, cosd(theta)] * 3.8 * 2
                end
            end
            residues[i] = AAResidue(aa_residue_name, [i])
        end
        new_chain = AAChain(""*chain_id_lib[1], rpad(mol_name, 6)[1:6], MOL_PROTEIN, [i for i = 1 : protein_length])

        chain0[1] = new_chain
    elseif args["strategy"] == "random-walk"
        println("Self-avoiding random walk not support yet.")
    end

    new_mol0 = AAMolecule(atom_names[1:protein_length], atom_coors[1:3, 1:protein_length], residues[1:protein_length], chain0)

    # ===============================
    # coarse graining from AAMolecule
    # ===============================
    force_field = ForceFieldCG(ff_pro, 1, 1, 0, 0, 0)
    if idr_model == 1
        args["modeling-options"] = Dict("IDR" => Dict("HPS_region" => "1 to $protein_length"))
    elseif idr_model == 2
        args["modeling-options"] = Dict("IDR" => Dict("KH_region" => "1 to $protein_length"))
    elseif idr_model == 3
        args["modeling-options"] = Dict("IDR" => Dict("AICG2p_IDR_local" => "1 to $protein_length",
                                                      "AICG2p_IDR_nonlocal" => "1 to $protein_length"))
    end
    args["cgconnect"] = true

    cg_top0, cg_conf0 = coarse_graining(new_mol0, force_field, args)

    write_grotop(cg_top0, mol_name, args)
    write_grocrd(cg_top0, cg_conf0, mol_name, args)
    write_pdb(cg_top0, cg_conf0, mol_name, args)
    write_psf(cg_top0, mol_name, args)

    return 0
end


function main()
    args = parse_commandline()
    make_cg_protein_structure(args)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
