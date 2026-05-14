#!/usr/bin/env julia

include("../../../src/lib/gcj.jl")

using Random, Dates
using Printf
using Statistics
using ArgParse

function complex_duplication(pdb_filename::String, n_chains::Int=20, n_trials::Int=100, d_min::Float64=2.0, i_try::Int=1, size_scale::Float64=1.0)

    # ======================
    # prepare data structure
    # ======================
    chains_final = []

    # initialize random number generator
    seed = Int(round(time()))
    rng = MersenneTwister(seed)

    # ======================================
    # read single chain coordinates from PDB
    # ======================================
    println("-------------------------------------------------------------")
    println("Reading single chain coordinates from $pdb_filename")
    pdb_lines = readlines(pdb_filename)
    n_atoms = 0
    for line in pdb_lines
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            n_atoms += 1
        end
    end
    single_chain_coors = zeros(Float64, (3, n_atoms))
    single_chain_line_head = []
    single_chain_line_tail = []
    i_atom = 0
    for line in pdb_lines
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            i_atom += 1
            # extract coordinates from the line
            x = parse(Float64, line[31:38])
            y = parse(Float64, line[39:46])
            z = parse(Float64, line[47:54])
            single_chain_coors[:, i_atom] = [x, y, z]
            line_head = line[1:30]    # store the line head
            line_tail = line[55:end]  # store the line tail
            push!(single_chain_line_head, line_head)  # store the line head
            push!(single_chain_line_tail, line_tail)  # store the line tail
        end
    end

    # move COM of single chain to origin
    println("Number of atoms in single chain: $n_atoms")
    println("Computing center of mass (COM) of the single chain...")
    com = mean(single_chain_coors, dims=2)
    single_chain_coors .-= com

    # -----------------------
    # determine system size L
    # -----------------------
    coor_1 = single_chain_coors[:, 1]
    coor_2 = single_chain_coors[:, end]
    chain_len_0 = compute_distance(coor_1, coor_2)

    L = chain_len_0 * size_scale


    # --------------------
    # make a box of size L
    # --------------------
    x_min = -L / 2
    x_max = L / 2
    y_min = -L / 2
    y_max = L / 2
    z_min = -L / 2
    z_max = L / 2

    # ------------------------------------------------------
    # divide the box into smaller cells for cell-linked list
    # ------------------------------------------------------
    println("Dividing the box into cells for efficient placement checks...")
    n_cells_x = Int(ceil(L / d_min / 5)) # 5 is an arbitrary factor to ensure enough cells
    n_cells_y = n_cells_x
    n_cells_z = n_cells_x

    # Create a 3D array (cell list) to hold atoms in each cell
    # cell_list = [Vector{Tuple{Int, Int}}() for i in 1:n_cells_x, j in 1:n_cells_y, k in 1:n_cells_z]
    # Use BitArray to minimize memory usage (each cell is 1 bit)
    cell_list = BitArray(undef, n_cells_x, n_cells_y, n_cells_z)
    fill!(cell_list, false)

    # Helper function to map coordinates to cell indices
    function get_cell_indices(x, y, z)
        ix = clamp(Int(floor((x - x_min) / (x_max - x_min) * n_cells_x)) + 1, 1, n_cells_x)
        iy = clamp(Int(floor((y - y_min) / (y_max - y_min) * n_cells_y)) + 1, 1, n_cells_y)
        iz = clamp(Int(floor((z - z_min) / (z_max - z_min) * n_cells_z)) + 1, 1, n_cells_z)
        return ix, iy, iz
    end

    # ==============================
    # put chains into simulation box
    # ==============================
    println("Placing $n_chains chains into the simulation box...")
    for i_chain in 1:n_chains
        println(" > chain $i_chain of $n_chains...")
        success = false
        for t in 1:n_trials
            println("   > trial $t of $n_trials...")
            # generate a random rotation matrix
            q = generate_random_rotation()

            # rotate the single chain coordinates
            coords_rot = q * single_chain_coors

            # generate a random translation vector within the box
            t_vec = randn(rng, 3) * (L / 10)  # normal distribution, stddev = L/10

            # translate the rotated coordinates
            coords_try = coords_rot .+ t_vec

            # ---------------------------------------------------------------------------------------------
            # check if the minimum distance between the new chain and existing chains is greater than d_min
            # ---------------------------------------------------------------------------------------------
            is_a_good_placement = true
            for atom_idx in 1:n_atoms
                x, y, z = coords_try[:, atom_idx]
                ix, iy, iz = get_cell_indices(x, y, z)
                if cell_list[ix, iy, iz]
                    is_a_good_placement = false
                    break
                end
            end

            if is_a_good_placement
                # Place all atoms from the new chain into cell_list (optional, if needed for overlap check)
                for atom_idx in 1:n_atoms
                    x, y, z = coords_try[:, atom_idx]
                    ix, iy, iz = get_cell_indices(x, y, z)
                    # push!(cell_list[ix, iy, iz], (i_chain, atom_idx))
                    cell_list[ix, iy, iz] = true
                end

                push!(chains_final, coords_try)
                success = true
                break
            end
        end

        if !success
            println("Failed to place chain $i_chain after $n_trials trials.")
        end
    end


    # ==================
    # output coordinates
    # ==================
    println("Outputting coordinates to files...")

    sysname = "__gen__" * pdb_filename[1:end-4] * "_20chains"
    # ---
    # PDB
    # ---
    newpdb_filename = @sprintf("%s_%02d.pdb", sysname, i_try)
    newpdb_file = open(newpdb_filename, "w")
    for i_chain in 1:n_chains
        for atom_idx in 1:n_atoms
            coor = chains_final[i_chain][:, atom_idx]
            line_head = single_chain_line_head[atom_idx]
            line_tail = single_chain_line_tail[atom_idx]
            chain_id = Char('A' + i_chain - 1)
            line_head_new = line_head[1:21] * chain_id * line_head[23:end]  # modify chain ID in the line head
            @printf(newpdb_file, "%s%8.3f%8.3f%8.3f%s\n", line_head_new, coor[1], coor[2], coor[3], line_tail)
        end
    end


    # ------
    # inpcrd
    # ------
    inpcrd_filename = @sprintf("%s_%02d.inpcrd", sysname, i_try)
    inpcrd_file = open(inpcrd_filename, "w")
    println(inpcrd_file, "20 chains generated from $pdb_filename")
    println(inpcrd_file, n_chains * n_atoms)
    for i_chain in 1:n_chains
        for atom_idx in 1:n_atoms
            coor = chains_final[i_chain][:, atom_idx]
            @printf(inpcrd_file, "%12.7f%12.7f%12.7f", coor[1], coor[2], coor[3])
            if atom_idx % 2 == 0
                @printf(inpcrd_file, "\n")
            end
        end
    end

end


# =============================
# Parsing Commandline Arguments
# =============================
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin

        "--top", "-t"
        help     = "Topology file name (gromacs style)."
        arg_type = String

        "--crd", "-c"
        help     = "Coordinate file name (gromacs style)."
        arg_type = String

        "--pdb", "-P"
        help     = "Coordinate file name (PDB format)."
        required = true
        arg_type = String

        "--output", "-o"
        help     = "Output file name."
        arg_type = String
        default  = "BIGMOL"

        "--num-chain", "-n"
        help     = "Number of copies."
        arg_type = Int
        default  = 1

        "--num-trial"
        help     = "Number of trials in randomization."
        arg_type = Int
        default  = 100

        "--num-output"
        help     = "Number of new structures."
        arg_type = Int
        default  = 1

        "--min-distance", "-d"
        help     = "Density (probability) of molecules."
        arg_type = Float64
        default  = 1.0

        "--size-scale", "-S"
        help     = "Box size scale factor."
        arg_type = Float64
        default  = 1.5

        "--debug"
        help = "DEBUG."
        action = :store_true
    end

    return parse_args(s)
end


if abspath(PROGRAM_FILE) == @__FILE__

    args = parse_commandline()


    top_filename = get(args, "top", "")
    crd_filename = get(args, "crd", "")
    pdb_filename = get(args, "pdb", "")

    sys_name     = get(args, "output", "new_mol")

    n_chains = get(args, "num-chain", 1)
    n_trials = get(args, "num-trial", 100)

    d_min = get(args, "min-distance", 2.0)
    s_sca = get(args, "size-scale", 2.0)

    n_try = get(args, "num-output", 1)

    for i_try in 1:n_try
        complex_duplication(pdb_filename, n_chains, n_trials, d_min, i_try, s_sca)
    end
end
