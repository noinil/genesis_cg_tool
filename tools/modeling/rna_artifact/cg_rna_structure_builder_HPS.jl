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
        help = "RNA sequence file."
        arg_type = String
        default = ""

        "--length"
        help = "Number of nucleotides in the RNA with random-sequence."
        arg_type = Int
        default = 100

        "--straightness"
        help = "Angle threshold to create a (non)-straight chain."
        arg_type = Float64
        default = 45.0

    end

    return parse_args(s)
end


function make_cg_RNA_structure(args)

    seq_name = get(args, "sequence", "")

    # non-straightness (because ideal straight chain could have problem...)
    threshold_angle = get( args, "straightness", 45.0)

    println("============================================================")

    AA_FULLNAME_DICT = Dict(
        'A' => "RA",
        'C' => "RC",
        'G' => "RG",
        'U' => "RU"
    )
    MASS_DICT = Dict(
        'A' => 329.200,
        'C' => 305.200,
        'G' => 345.200,
        'U' => 306.200
    )

    # ========================================================
    # RNA sequence (read from file or generate random one)
    # ========================================================
    if length(seq_name) > 0
        println("> Open sequence file:", seq_name)

        mol_name = split(basename( seq_name ), '.')[1] * "_cg"

        RNA_seqence = ""
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
                error("Only support single-chain RNA!")
            end
            seq = strip(line)
            if length(seq) == 0
                continue
            end
            for b in seq
                if ! haskey(AA_FULLNAME_DICT, b)
                    error("Wrong RNA sequence!")
                end
            end
            RNA_seqence *= seq
        end
        RNA_length = length(RNA_seqence)
    else
        println("> Generating random RNA sequence:")

        mol_name = "random_RNA_cg"

        RNA_length = args["length"]
        RNA_seqence = randstring("ACGU", RNA_length)
    end
    println("> RNA sequence:   ( Length: $RNA_length ) ")
    println("> ", RNA_seqence)

    # ===================
    # Preparing RNA chain
    # ===================
    atom_types = Vector{String}(undef, RNA_length)
    atom_names = Vector{String}(undef, RNA_length)
    resi_names = Vector{String}(undef, RNA_length)
    atom_masss = zeros(RNA_length)
    atom_coors = zeros(3, RNA_length)

    theta = rand() * threshold_angle
    for i in 1 : RNA_length
        short_name = RNA_seqence[i]
        residue_name = AA_FULLNAME_DICT[short_name]
        a_mass = MASS_DICT[short_name]

        atom_types[i] = residue_name
        atom_names[i] = "RP"
        resi_names[i] = residue_name
        atom_masss[i] = a_mass

        if i == 1
            atom_coors[:, i] = [0, 0, 0]
        else
            if mod(i, 2) == 0
                theta = rand() * threshold_angle
                phi   = rand() * 360
                atom_coors[:, i] = atom_coors[:, i - 1] + [sind(theta) * cosd(phi), sind(theta) * sind(phi), cosd(theta)] * 5.0
            else
                atom_coors[:, i] = atom_coors[:, i - 2] + [0, 0, cosd(theta)] * 5.0 * 2
            end
        end
    end

    # ===============
    # Output topology
    # ===============
    # ---
    # top
    # ---
    top_fname = mol_name * ".top"
    top_file = open(top_fname, "w")

    @printf(top_file, "#include \"./param/atom_types.itp\" \n")
    @printf(top_file, "#include \"./param/flexible_local_angle.itp\" \n")
    @printf(top_file, "#include \"./param/flexible_local_dihedral.itp\" \n")
    @printf(top_file, "#include \"./param/pair_energy_MJ_96.itp\" \n\n")
    @printf(top_file, "#include \"./itp/%s.itp\" \n\n", mol_name)

    @printf(top_file, "[ system ] \n")
    @printf(top_file, "%s \n\n", mol_name)

    @printf(top_file, "[ molecules ] \n")
    @printf(top_file, "%s 1 \n\n", mol_name)

    @printf(top_file, "; [ cg_ele_chain_pairs ]\n")
    @printf(top_file, "; ON 1 - 2 : 3 - 4\n")

    close(top_file)

    # ---
    # itp
    # ---
    itp_fname = mol_name * ".itp"
    itp_file = open(itp_fname, "w")

    @printf(itp_file, "[ moleculetype ]\n")
    @printf(itp_file, "%s   3 \n", mol_name)

    @printf(itp_file, "\n[ atoms ]\n")
    @printf(itp_file, "; +INFO+ CHAIN:      1     SEGNAME: RAND_RNA \n")
    for i in 1 : RNA_length
        @printf(itp_file, "%10d%5s%10d%5s%5s%5d %8.3f %8.3f\n",
                i, atom_types[i], i, resi_names[i], atom_names[i], 1, -1.000, atom_masss[i])
    end

    @printf(itp_file, "\n[ bonds ]\n")
    for i in 1 : RNA_length - 1
        @printf(itp_file, "%10d%10d%5d%18.4E%18.4E\n",
                i, i + 1, 1, 0.5, 2000.0)
    end

    @printf(itp_file, "\n[ cg_IDR_HPS_region ]\n")
    @printf(itp_file, "%10d %10d\n", 1, RNA_length)

    close(itp_file)

    # ==================
    # Output coordinates
    # ==================
    crd_fname = mol_name * ".gro"
    crd_file = open(crd_fname, "w")

    @printf(crd_file, "CG model for %s, t = 0.000 \n", mol_name)
    @printf(crd_file, " %10d \n", RNA_length)

    for i in 1 : RNA_length
        @printf(crd_file, "%5d%5s%5s%5d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n",
                i, resi_names[i], atom_names[i], i,
                atom_coors[1, i] * 0.1,
                atom_coors[2, i] * 0.1,
                atom_coors[3, i] * 0.1,
                0.0, 0.0, 0.0)
    end
    @printf(crd_file, "%15.4f%15.4f%15.4f \n\n", 0.0, 0.0, 0.0)

    close(crd_file)

    return 0
end


function main()
    args = parse_commandline()
    make_cg_RNA_structure(args)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
