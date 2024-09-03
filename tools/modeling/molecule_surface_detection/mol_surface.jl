#!/usr/bin/env julia

using Random
using Printf
using ArgParse

include("../../../src/lib/gcj.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin

        "--GROCRD", "-c"
        help = "Protein structure file in the format of GROCRD."
        arg_type = String
        default = ""

        "--PDB", "-s"
        help = "Protein structure file in the format of PDB."
        arg_type = String
        default = ""

        "--mmCIF"
        help = "Use mmCIF format PDB file as input."
        action = :store_true

        "--tip-size", "-t"
        help = "Tip radius (Å)."
        arg_type = Float64
        default = 5.0

        "--grid-size", "-g"
        help = "Grid size (Å)."
        arg_type = Float64
        default = 5.0

        "--padding", "-b"
        help = "Padding distance (Å)."
        arg_type = Float64
        default = 15.0

        "--verbose", "-v"
        help = "Output more information."
        action = :store_true

    end

    return parse_args(s)
end


function determine_mol_surface(atom_coors::Array{<:Real,2}, atom_radii::Vector{<:Real}, tip_radius::Float64, grid_size::Float64, padding::Float64)

    n_atoms = size(atom_coors)[2]

    # ===============
    # box preparation
    # ===============
    # ------------------------
    # find min-max of molecule
    # ------------------------
    mol_coor_min     = minimum(atom_coors, dims=2)[1:3]
    mol_coor_max     = maximum(atom_coors, dims=2)[1:3]
    # geometric center
    mol_coor_cent    = 0.5 .* (mol_coor_max .+ mol_coor_min)
    # molecule size
    mol_size         = mol_coor_max .- mol_coor_min
    # move molecule to box center (origin)
    atom_coors_inbox = atom_coors .- mol_coor_cent
    # output
    println("> Molecule min:", mol_coor_min)
    println("> Molecule max:", mol_coor_max)

    # ------------------
    # determine box size
    # ------------------
    # get the maximum of atom radius
    atom_radii_max = maximum(atom_radii)
    # set a safe padding
    if padding < 2 * tip_radius + atom_radii_max
        padding = 2 * tip_radius + atom_radii_max
        @printf(" Padding distance was too small... Reset to %8.3f \n", padding)
    end
    box_size = mol_size .+ 2 * padding
    box_geo_min = box_size .* (-0.5)
    box_geo_max = box_size .* 0.5
    # output info
    println("> Box size:", box_size)

    # =====================
    # divide box into cells
    # =====================
    # estimate number of cells and cell size
    n_cell = Int.(ceil.(box_size ./ grid_size))
    cell_size = box_size ./ n_cell
    # output info
    println("> Grid number:", n_cell)
    println("> Grid size:", cell_size)

    # ================
    # find empty cells
    # ================
    cell_occupancy         = zeros(Int8, prod(n_cell))    # stat of cell: 0-empty, 1-occupied
    cell_closest_atom      = zeros(Int, prod(n_cell))     # idx of the closest atom
    cell_closest_distance  = zeros(Float64, prod(n_cell)) # dist from the closest atom
    cell_closest_distance .= 1e20                         # initialized as 1e20

    # convert 3d index to 1d index
    function cell_idx_3_to_1(idx::Vector{Int})
        return (idx[3] - 1) * n_cell[2] * n_cell[1] + (idx[2] - 1) * n_cell[1] + idx[1]
    end

    # convert 1d index to 3d index
    function cell_idx_1_to_3(idx::Int)
        iyz, ix = divrem(idx - 1, n_cell[1])
        iz,  iy = divrem(iyz,     n_cell[2])
        return [ix + 1, iy + 1, iz + 1]
    end

    for i_atom in 1:n_atoms
        r_atom = atom_radii[i_atom]

        # cell of atom i
        x_atom    = atom_coors_inbox[:, i_atom]
        i_cell    = Int.(ceil.((x_atom .- box_geo_min) ./ cell_size))
        i_cell_1d = cell_idx_3_to_1(i_cell)
        cell_occupancy[i_cell_1d] = 1

        # -------------------------------
        # set occupancy of neighbor cells
        # -------------------------------
        m_i = Int.(round.((r_atom + tip_radius) ./ cell_size))
        for dx in -m_i[1] : m_i[1]
            for dy in -m_i[2] : m_i[2]
                for dz in -m_i[3] : m_i[3]
                    j_cell    = i_cell .+ [dx, dy, dz]
                    j_cell_1d = cell_idx_3_to_1(j_cell)
                    # if j_cell is out of range...
                    if any(<(1), j_cell) || any(j_cell > n_cell)
                        continue
                    end
                    if cell_occupancy[j_cell_1d] == 1
                        continue
                    end
                    # coordinte of neighboring cell center
                    j_cell_center = box_geo_min + (j_cell .- 0.5) .* cell_size
                    # calculate distance
                    d_ij = compute_distance(x_atom, j_cell_center)
                    # compare with r_tip + r_atom
                    if d_ij <= r_atom + tip_radius
                        cell_occupancy[j_cell_1d] = 1
                    end
                    # compare the distance
                    if d_ij < cell_closest_distance[j_cell_1d]
                        cell_closest_distance[j_cell_1d] = d_ij
                        cell_closest_atom[j_cell_1d]     = i_atom
                    end
                end
            end
        end
        # end loop over neighbor cells
    end
    # end setting occupancy


    # =======================
    # determine surface atoms
    # =======================
    atom_visibility = zeros(Int8, n_atoms)
    for i_cell_1d in findall(<(1), cell_occupancy)
        i_atom = cell_closest_atom[i_cell_1d]
        d_atom = cell_closest_distance[i_cell_1d]
        if i_atom > 0 && d_atom < 1e10
            atom_visibility[i_atom] = 1
        end
    end

    return atom_visibility
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = parse_commandline()

    pdb_name = get(args, "PDB", "")
    crd_name = get(args, "GROCRD", "")
    is_mmCIF = get(args, "mmCIF", false)
    verbose  = get(args, "verbose", false)
    tsize    = args["tip-size"]
    gsize    = args["grid-size"]
    psize    = args["padding"]

    if length(pdb_name) > 0
        if is_mmCIF
            cif_data = read_mmCIF(pdb_name)
            aa_molecule = mmCIF_to_AAMolecule(cif_data)
        else
            aa_molecule = read_PDB(pdb_name)
        end
        aa_num_atom    = length(aa_molecule.atom_names)
        aa_num_residue = length(aa_molecule.residues)
        aa_num_chain   = length(aa_molecule.chains)

        if verbose
            println("          > Number of atoms    : $(aa_num_atom)")
            println("          > Number of residues : $(aa_num_residue)")
            println("          > Number of chains   : $(aa_num_chain)")
        end

        # --------------------
        # preparing atom radii
        # --------------------
        atom_radii = zeros(Float64, aa_num_atom)
        atom_radius_lib = Dict(
            'H'=>1.2,
            'C'=>1.7,
            'O'=>1.52,
            'N'=>1.55,
            'P'=>1.8,
            'S'=>1.8
        )
        for j_atom in 1:aa_num_atom
            a_name = aa_molecule.atom_names[j_atom]
            r_atom = atom_radius_lib[a_name[1]]
            atom_radii[j_atom] = r_atom
        end

        # --------
        # go go go
        # --------
        surface_atoms = determine_mol_surface(aa_molecule.atom_coors, atom_radii, tsize, gsize, psize)



        # --------------------------
        # post-surface-determination
        # --------------------------
        # atomic surface and residual surface
        surf_res_list = zeros(Int8, aa_num_residue)
        for j_res in 1:aa_num_residue
            for j_atom in aa_molecule.residues[j_res].atoms
                if surface_atoms[j_atom] > 0
                    surf_res_list[j_res] = 1
                    break
                end
            end
            if surf_res_list[j_res] < 1
                continue
            end
            for j_atom in aa_molecule.residues[j_res].atoms
                surface_atoms[j_atom] = 1
            end
        end

        # test output
        new_pdb_name = @sprintf("__test_%s_%3.1f_%3.1f_%3.1f_.pdb", pdb_name[1:end-4], tsize, gsize, psize)
        new_pdb = open(new_pdb_name, "w")
        # re-output new pdb
        k_atom = 0
        for original_pdb_line in readlines(pdb_name)
            if startswith(original_pdb_line, "ATOM  ")
                global k_atom += 1
                occup = surface_atoms[k_atom]
                @printf(new_pdb, "%s%8.3f%s\n", original_pdb_line[1:54], occup, original_pdb_line[61:end])
            end
        end
        close(new_pdb)
    end
end
