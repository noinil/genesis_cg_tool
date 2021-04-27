###############################################################################
#                            ____   ____   ____                               #
#                           |  _ \ |  _ \ | __ )                              #
#                           | |_) || | | ||  _ \                              #
#                           |  __/ | |_| || |_) |                             #
#                           |_|    |____/ |____/                              #
#                                                                             #
###############################################################################

# =============
# Function list
#
# function parse_PDB_line(pdb_line::AbstractString)
# function write_PDB_line(io::IO, aline::PDBLine)
# function read_PDB(pdb_name::AbstractString)
# function write_pdb(top::GenTopology, conf::Conformation, system_name::AbstractString, args::Dict{String, <:Any}=Dict{String, Any}())
# =============

using Printf

struct PDBLine
    atom_serial::Int     # line[7:11]
    atom_name::String    # line[13:16]
    residue_name::String # line[18:21]
    chain_id::String     # line[22]
    residue_serial::Int  # line[23:26]
    coor_x::Float64      # line[31:38]
    coor_y::Float64      # line[39:46]
    coor_z::Float64      # line[47:54]
    occupancy::Float64   # line[55:60]
    tempfactor::Float64  # line[61:66]
    segment_id::String   # line[67:76]
    element_name::String # line[77:78]
    charge::Float64      # line[79:80]
end

function parse_PDB_line(pdb_line::AbstractString)
    atom_serial     = 0
    try
        atom_serial = parse(Int, pdb_line[7:11])
    catch
        atom_serial = 0
        # println("WARNING! Error in reading atom serial in PDB!")
    end
    atom_name       = strip(pdb_line[13:16])
    residue_name    = strip(pdb_line[18:21])
    chain_id        = pdb_line[22:22]
    residue_serial  = parse(Int, pdb_line[23:26])
    coor_x          = parse(Float64, pdb_line[31:38])
    coor_y          = parse(Float64, pdb_line[39:46])
    coor_z          = parse(Float64, pdb_line[47:54])
    occupancy       = 0.0
    try
        occupancy   = parse(Float64, pdb_line[55:60])
    catch
        occupancy   = 0.0
    end
    tempfactor      = 0.0
    try
        tempfactor  = parse(Float64, pdb_line[61:66])
    catch
        tempfactor  = 0.0
    end
    segment_id      = strip(pdb_line[67:76])
    element_name    = strip(pdb_line[77:78])
    charge          = 0.0
    try
        charge      = parse(Float64, pdb_line[79:80])
    catch
        charge      = 0.0
    end
    new_pdb_data = PDBLine(atom_serial, atom_name,
                           residue_name, chain_id, residue_serial,
                           coor_x, coor_y, coor_z,
                           occupancy, tempfactor, segment_id,
                           element_name, charge)
    return new_pdb_data
end

function write_PDB_line(io::IO, aline::PDBLine)
    @printf(io, "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s%2d \n",
            aline.atom_serial,
            aline.atom_name,
            rpad(aline.residue_name, 4),
            aline.chain_id,
            aline.residue_serial,
            aline.coor_x,
            aline.coor_y,
            aline.coor_z,
            aline.occupancy,
            aline.tempfactor,
            aline.segment_id,
            aline.element_name,
            Int(round(aline.charge)))
end



# ===========
# Read PDB!!!
# ===========

function read_PDB(pdb_name::AbstractString)
    aa_pdb_lines = []

    # =================================
    # Step 1: Determine number of atoms
    # =================================

    aa_num_atom = 0
    for line in eachline(pdb_name)
        if startswith(line, "ATOM")
            push!(aa_pdb_lines, rpad(line, 80))
            aa_num_atom += 1
        elseif startswith(line, "TER") || startswith(line, "END")
            push!(aa_pdb_lines, line)
        end
    end

    # ==========================
    # Data structures for output
    # ==========================

    aa_atom_name  = fill("    ",       aa_num_atom)
    aa_coor       = zeros(Float64, (3, aa_num_atom))
    aa_residues   = []
    aa_chains     = []

    # ---------------
    # Local variables
    # ---------------

    i_atom        = 0
    i_resid       = 0
    curr_resid    = NaN
    curr_chain    = NaN
    curr_rname    = "    "
    residue_name  = "    "
    chain_id      = "?"
    tmp_res_atoms = []
    tmp_chain_res = []
    segment_id    = " "

    # ========================================
    # Step 2: Add atoms to residues and chains
    # ========================================

    for line in aa_pdb_lines
        if startswith(line, "TER") || startswith(line, "END")
            if length(tmp_res_atoms) > 0
                push!(aa_residues, AAResidue(residue_name, tmp_res_atoms))
                tmp_res_atoms = []
            end
            if length(tmp_chain_res) > 0
                # -----------------------------
                # Determine chain molecule type
                # -----------------------------
                mol_type = -1
                for i_res in tmp_chain_res
                    res_name = aa_residues[i_res].name
                    tmp_mol_type = MOL_OTHER
                    if in(res_name, RES_NAME_LIST_PROTEIN)
                        tmp_mol_type = MOL_PROTEIN
                    elseif in(res_name, RES_NAME_LIST_DNA)
                        tmp_mol_type = MOL_DNA
                    elseif in(res_name, RES_NAME_LIST_RNA)
                        tmp_mol_type = MOL_RNA
                    elseif haskey(RES_NAME_RNA_DICT, res_name) || haskey(RES_NAME_DNA_DICT, res_name)
                        tmp_mol_type = MOL_DNA
                        for i_atom in aa_residues[i_res].atoms
                            atom_name = aa_atom_name[i_atom]
                            if atom_name == "O2'"
                                tmp_mol_type = MOL_RNA
                                break
                            end
                        end
                    end
                    if mol_type == -1
                        mol_type = tmp_mol_type
                    elseif tmp_mol_type != mol_type
                        errmsg = @sprintf("BUG: Inconsistent residue types in chain ID - %s residue - %d : %s ",
                                          chain_id,
                                          i_res,
                                          res_name)
                        error(errmsg)
                    end
                end
                # --------------------------------------
                # chain mol type determination ends here
                # --------------------------------------

                push!(aa_chains, AAChain(chain_id, segment_id, mol_type, tmp_chain_res))
                tmp_chain_res = []
            end
            continue
        end

        new_pdb_data = parse_PDB_line(line)

        i_atom        += 1
        atom_name      = new_pdb_data.atom_name
        residue_name   = new_pdb_data.residue_name
        chain_id       = new_pdb_data.chain_id
        residue_serial = new_pdb_data.residue_serial
        coor_x         = new_pdb_data.coor_x
        coor_y         = new_pdb_data.coor_y
        coor_z         = new_pdb_data.coor_z
        segment_id     = new_pdb_data.segment_id

        aa_atom_name[i_atom] = atom_name
        aa_coor[1, i_atom]   = coor_x
        aa_coor[2, i_atom]   = coor_y
        aa_coor[3, i_atom]   = coor_z

        if residue_serial != curr_resid
            i_resid += 1
            push!(tmp_chain_res, i_resid)
            curr_resid = residue_serial
            if length(tmp_res_atoms) > 0
                push!(aa_residues, AAResidue(curr_rname, tmp_res_atoms))
                tmp_res_atoms = []
            end
            curr_rname = residue_name
        end

        push!(tmp_res_atoms, i_atom)
    end

    new_molecule = AAMolecule(aa_atom_name, aa_coor, aa_residues, aa_chains)

    return new_molecule

end


# =============
# Output CG PDB
# =============
function write_pdb(top::GenTopology, conf::Conformation, system_name::AbstractString, args::Dict{String, <:Any}=Dict{String, Any}())

    verbose = get(args, "verbose", false)

    pdb_name        = system_name * ".pdb"
    pdb_file        = open(pdb_name, "w")

    do_output_cgconnect = get(args, "cgconnect", false)

    num_particles  = conf.num_particle
    is_huge_system = num_particles > 9999

    chain_id_set   = "_ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"
    tmp_chain_id   = 0
    tmp_seg_name   = ""
    real_chain_id  = 1
    for i_bead in 1 : num_particles
        i_chain = top.top_atoms[i_bead].chain_id
        i_sname = top.top_atoms[i_bead].seg_name
        if i_chain != tmp_chain_id || i_sname != tmp_seg_name
            if tmp_chain_id > 0
                print(pdb_file, "TER\n")
                real_chain_id += 1
            end
            tmp_chain_id = i_chain
            tmp_seg_name = i_sname
        end
        resid_index_tmp = top.top_atoms[i_bead].residue_index

        @printf(pdb_file,
                "ATOM  %5d %4s%1s%4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s%2s \n",
                i_bead % 100000,
                top.top_atoms[i_bead].atom_name,
                ' ',
                rpad( top.top_atoms[i_bead].residue_name, 4 ),
                chain_id_set[mod(real_chain_id, 63) + 1],
                resid_index_tmp % 10000,
                ' ',
                conf.coors[1 , i_bead],
                conf.coors[2 , i_bead],
                conf.coors[3 , i_bead],
                0.0,
                0.0,
                top.top_atoms[i_bead].seg_name,
                "",
                "")
    end
    print(pdb_file,"TER\n")

    cg_pdb_cnct_line = "CONECT%5d%5d \n"
    if do_output_cgconnect
        for bond in top.top_bonds
            @printf(pdb_file, "CONECT%5d%5d \n", bond.i, bond.j)
        end
    end

    print(pdb_file,"END\n")
    print(pdb_file,"\n")

    close(pdb_file)

    if verbose
        println(">           ... .pdb (CG) : DONE!")
    end

end


