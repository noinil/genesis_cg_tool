#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/biomath.jl")
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/constants.jl")
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/topology.jl")
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/conformation.jl")
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/parser_crd.jl")
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/parser_top.jl")

function main()

    NUM_X = 3
    NUM_Y = 3
    NUM_Z = 3

	CELL_SIZE_X = 100
	CELL_SIZE_Y = 100
	CELL_SIZE_Z = 100

    mol_name = "WHATEVER"

    mol_top = read_grotop("./mol_cg.top")
    mol_crd = read_grocrd("./gro/mol_cg.gro")

    println("System name:", mol_top.system_name)
    println("Number of particles in top:", mol_top.num_atom)

    # ==============================
    # move single molecule to origin
    # ==============================
    mol_center = centroid(mol_crd.coors)
    mol_orig_coors = mol_crd.coors .- mol_center

    # ==================
    # make copies of mol
    # ==================
    num_copies = NUM_X * NUM_Y * NUM_Z
    total_num_particles = mol_top.num_atom * num_copies

    gro_name = @sprintf("%s_mul_%d_%d_%d.gro", mol_name, NUM_X, NUM_Y, NUM_Z)
    gro_file = open(gro_name, "w")

    @printf(gro_file, "CG model %s, nx: %d, ny: %d, nz: %d, t = %16.3f \n", mol_top.system_name, NUM_X, NUM_Y, NUM_Z, 0)
    @printf(gro_file, "%12d \n", total_num_particles)

    @printf("Duplicated system has %d x %d x %d = %d copies, in total %d atoms \n", NUM_X, NUM_Y, NUM_Z, num_copies, total_num_particles)

    i_bead_global = 0
    for ix in 1:NUM_X
        for iy in 1:NUM_Y
            for iz in 1:NUM_Z
                shift_x = (ix - 1) * CELL_SIZE_X
                shift_y = (iy - 1) * CELL_SIZE_Y
                shift_z = (iz - 1) * CELL_SIZE_Z
                for i_bead in 1 : mol_top.num_atom
                    i_bead_global += 1
                    @printf(gro_file, "%5d%5s%5s%5d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n",
                            mol_top.top_atoms[i_bead].residue_index,
                            mol_top.top_atoms[i_bead].residue_name,
                            mol_top.top_atoms[i_bead].atom_name,
                            i_bead_global % 100000,
                            ( mol_orig_coors[1,i_bead] + shift_x ) * 0.1,
                            ( mol_orig_coors[2,i_bead] + shift_y ) * 0.1,
                            ( mol_orig_coors[3,i_bead] + shift_z ) * 0.1,
                            0.0, 0.0, 0.0)
                end
            end
        end
    end
    @printf(gro_file, "%15.4f%15.4f%15.4f \n\n", 0.0, 0.0, 0.0)

    close(gro_file)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
