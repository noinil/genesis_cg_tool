#!/usr/bin/env julia

using Printf
using ArgParse

include("../../src/lib/biomath.jl")
include("../../src/lib/topology.jl")
include("../../src/lib/conformation.jl")
include("../../src/lib/parser_top.jl")
include("../../src/lib/parser_crd.jl")

function main(args)

    top_filename = get(args, "top", "")
    crd_filename = get(args, "crd", "")

    # ================================
    # Read in topology and coordinates
    # ================================

    mytop = read_grotop(top_filename)
    mycrd = read_grocrd(crd_filename)

    # ===============
    # Structure check
    # ===============

    # -----------
    # bond length
    # -----------
    println("====================================================================================================")
    for bond in mytop.top_bonds
        if bond.r0 < 0.3
            @printf("Short bond length (< 3.0Å): %10d %3s %10d %3s - %6.3fÅ \n",
                    bond.i, mytop.top_atoms[bond.i].atom_type,
                    bond.j, mytop.top_atoms[bond.j].atom_type,
                    bond.r0 * 10)
        end
    end

    # -----
    # angle
    # -----
    println("====================================================================================================")
    for angle in mytop.top_angles
        if angle.function_type < 20 && angle.a0 > 165.0
            @printf("Large angle (close to π): %10d %3s %10d %3s %10d %3s - %6.3f° \n",
                    angle.i, mytop.top_atoms[angle.i].atom_type,
                    angle.j, mytop.top_atoms[angle.j].atom_type,
                    angle.k, mytop.top_atoms[angle.k].atom_type,
                    angle.a0)
        elseif angle.function_type == 21 && angle.a0 > 0.75
            coor_i = mycrd.coors[:, angle.i]
            coor_j = mycrd.coors[:, angle.j]
            coor_k = mycrd.coors[:, angle.k]
            ang = compute_angle(coor_i, coor_j, coor_k)
            @printf("Large angle (close to π): %10d %3s %10d %3s %10d %3s - %6.3f° \n",
                    angle.i, mytop.top_atoms[angle.i].atom_type,
                    angle.j, mytop.top_atoms[angle.j].atom_type,
                    angle.k, mytop.top_atoms[angle.k].atom_type,
                    ang)
        end
    end

    # ---------------
    # native contacts
    # ---------------
    println("====================================================================================================")
    for contact in mytop.top_pairs
        if contact.r0 < 0.4
            @printf("Short contact (< 4.0Å): %10d %3s %10d %3s - %6.3fÅ \n",
                    contact.i, mytop.top_atoms[contact.i].atom_type,
                    contact.j, mytop.top_atoms[contact.j].atom_type,
                    contact.r0 * 10)
        end
    end

    println("====================================================================================================")
end

# =============================
# Parsing Commandline Arguments
# =============================
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin

        "--top", "-t"
        help     = "Topology file name (gromacs style)."
        required = true
        arg_type = String

        "--crd", "-c"
        help     = "Coordinate file name (gromacs style)."
        required = true
        arg_type = String

        "--debug"
        help = "DEBUG."
        action = :store_true
    end

    return parse_args(s)
end



if abspath(PROGRAM_FILE) == @__FILE__

    args = parse_commandline()

    main(args)

end
