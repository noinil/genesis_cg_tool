#!/usr/bin/env julia

using Random
using Printf
using ArgParse

include("../../../src/lib/biomath.jl")
include("../../../src/lib/molecule.jl")
include("../../../src/lib/topology.jl")
include("../../../src/lib/constants.jl")
include("../../../src/lib/selection.jl")
include("../../../src/lib/coarse_graining_subroutines.jl")
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

        "--copy"
        help = "Number of chains created."
        arg_type = Int
        default = 1

        "--dx"
        help = "Delta x of the chain grid."
        arg_type = Float64
        default = 10.0

        "--dy"
        help = "Delta y of the chain grid."
        arg_type = Float64
        default = 10.0

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

    grid_dx = get(args, "dx", 10.0)
    grid_dy = get(args, "dy", 10.0)
    mol_num = get(args, "copy", 1)
    mol_grid_size = floor(Int, sqrt(mol_num))

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
    atom_names = Vector{String}(undef, protein_length * mol_num)
    atom_coors = zeros(3, protein_length * mol_num)
    residues   = Vector{AAResidue}(undef, protein_length * mol_num)
    chains     = Vector{AAChain}(undef, mol_num)
    chain0     = Vector{AAChain}(undef, 1)

    if args["strategy"] == "straight"
        for l in 1 : mol_num
            i_offset = (l - 1) * protein_length
            theta = rand() * threshold_angle
            for i in 1 : protein_length
                j = i + i_offset

                aa_short_name = protein_seqence[i]
                aa_residue_name = AA_FULLNAME_DICT[aa_short_name]
                
                atom_names[j] = "CA"
                if i == 1
                    (m, n) = fldmod(l - 1, mol_grid_size)
                    atom_coors[:, j] = [m * grid_dx, n * grid_dy, 0]
                else
                    if mod(i, 2) == 0
                        theta = rand() * threshold_angle
                        phi   = rand() * 360
                        atom_coors[:, j] = atom_coors[:, j - 1] + [sind(theta) * cosd(phi), sind(theta) * sind(phi), cosd(theta)] * 3.8
                    else
                        atom_coors[:, j] = atom_coors[:, j - 2] + [0, 0, cosd(theta)] * 3.8 * 2
                    end
                end
                residues[j] = AAResidue(aa_residue_name, [j])
            end
            new_chain = AAChain(chain_id_lib[mod(l, 63) + 1], rpad(mol_name, 6)[1:6], MOL_PROTEIN, [i + i_offset for i = 1 : protein_length])
            chains[l] = new_chain
        end
        chain0[1] = chains[1]
    elseif args["strategy"] == "random-walk"
        println("Self-avoiding random walk not support yet.")
    end

    new_mol0 = AAMolecule(atom_names[1:protein_length], atom_coors[1:3, 1:protein_length], residues[1:protein_length], chain0)
    new_mols = AAMolecule(atom_names, atom_coors, residues, chains)

    # ===============================
    # coarse graining from AAMolecule
    # ===============================
    force_field = ForceFieldCG(ff_pro, 1, 1, 0, 0, 0)
    args["modeling-options"] = Dict("IDR" => Dict("HPS_region" => "1 to $protein_length"))
    args["cgconnect"] = true

    cg_top0, cg_conf0 = coarse_graining(new_mol0, force_field, args)
    cg_tops, cg_confs = coarse_graining(new_mols, force_field, args)

    write_grotop(cg_top0, mol_name, args)
    write_grocrd(cg_tops, cg_confs, mol_name, args)
    write_pdb(cg_tops, cg_confs, mol_name, args)
    write_psf(cg_tops, mol_name, args)

    return 0
end


function main()
    args = parse_commandline()
    make_cg_protein_structure(args)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
