###############################################################################
#      ____  ____  ____          __                    ____ ___ _____         #
#     |  _ \|  _ \| __ )_  __   / / __ ___  _ __ ___  / ___|_ _|  ___|        #
#     | |_) | | | |  _ \ \/ /  / / '_ ` _ \| '_ ` _ \| |    | || |_           #
#     |  __/| |_| | |_) >  <  / /| | | | | | | | | | | |___ | ||  _|          #
#     |_|   |____/|____/_/\_\/_/ |_| |_| |_|_| |_| |_|\____|___|_|            #
###############################################################################

###############################################################################
# Function list
#
# read_mmCIF(cif_name::AbstractString, args::Dict{String, <:Any}=Dict{String, Any}())
# write_mmCIF(top::GenTopology, conf::Conformation, system_name::AbstractString, args::Dict{String, <:Any}=Dict{String, Any}())
###############################################################################

using Printf

struct mmCIF_Atom_Site_Info     # NOT REALLY USED
    group_PDB::String           # "ATOM" or "HETATM"
    id::Int                     # atom serial
    type_symbol::String         # atom type
    label_atom_id::String       # atom name
    label_alt_id::String        # identifier for alternative site
    label_comp_id::String       # 3-char residue name
    label_asym_id::String       # chain or asymmetric identifier
    label_entity_id::Int        # chain identifier
    label_seq_id::Int           # residue number
    pdbx_PDB_ins_code::String   # PDB insertion code
    Cartn_x::Float64            # coordinate x
    Cartn_y::Float64            # coordinate y
    Cartn_z::Float64            # coordinate z
    occupancy::Float64          # occupancy
    B_iso_or_equiv::Float64     # B factor
    pdbx_formal_charge::Float64 # charge
    auth_seq_id::Int            # AUTHOR residue number
    auth_comp_id::String        # AUTHOR 3-char residue name
    auth_asym_id::String        # AUTHOR chain identifier
    auth_atom_id::String        # AUTHOR atom name
    pdbx_PDB_model_num::Int     # model number
end

# ===========
# Read mmCIF!
# ===========
# TODO: multiline CHARACTER values
function mmcif_split(s::AbstractString)
    new_words = []

    word_tmp = ""
    is_substring = false
    for c in s
        if c == '\'' || c == '\"'
            is_substring = !is_substring
            continue
        end
        if is_substring
            word_tmp *= c
        else
            if c == ' '
                if length(word_tmp) > 0
                    push!(new_words, word_tmp)
                end
                word_tmp = ""
            else
                word_tmp *= c
            end
        end
    end
    if length(word_tmp) > 0
        push!(new_words, word_tmp)
    end
    return new_words[:]
end
function read_mmCIF(cif_name::AbstractString, args::Dict{String, <:Any}=Dict{String, Any}())

    # ==================================
    # The only structure is a dictionary
    # ==================================
    data_attributes = Dict()

    # ============
    # arguments...
    # ============
    verbose = get(args, "verbose", false)

    # ==================
    # local variables...
    # ==================
    TABULAR_FLAG = false
    TABULAR_KEYS = []
    MULTILINE_CHAR_FLAG = false

    # ===========
    # let's go...
    # ===========
    for line in eachline(cif_name)
        if startswith(line, "data")
            data_attributes["entry_name"] = strip(line)
            continue
        end
        if startswith(line, '#')
            continue
        end
        if startswith(line, ';')
            MULTILINE_CHAR_FLAG = ! MULTILINE_CHAR_FLAG
        end
        if MULTILINE_CHAR_FLAG
            continue
        end
        if strip(line) == "loop_"
            TABULAR_FLAG = true
            TABULAR_KEYS = []
            continue
        end
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            words = split(line)
        else
            words = mmcif_split(line)
        end
        # ---------
        # key-value
        # ---------
        if words[1][1] == '_' && length(words) > 1
            cif_key = words[1]
            cif_value = words[2]
            data_attributes[cif_key] = cif_value
            continue
        end
        # ------------
        # tabular keys
        # ------------
        if words[1][1] == '_' && length(words) == 1 && TABULAR_FLAG
            cif_key_tabu = words[1]
            push!(TABULAR_KEYS, cif_key_tabu)
            data_attributes[cif_key_tabu] = []
            continue
        end
        # --------------
        # tabular values
        # --------------
        if words[1][1] != '_'
            TABULAR_FLAG = false
            for ( i, w ) in enumerate(words)
                push!(data_attributes[TABULAR_KEYS[i]], w)
            end
        end
    end

    return data_attributes
end


function mmCIF_to_AAMolecule(cif_data::Dict)
    # ==============================================
    # Extract information frmo the _atom_site block!
    # ==============================================
    aa_num_atom = count(flag->(flag=="ATOM"), cif_data["_atom_site.group_PDB"])

    # ==============================
    # Data structures for AAMolecule
    # ==============================
    aa_atom_name  = fill("    ",       aa_num_atom)
    aa_coor       = zeros(Float64, (3, aa_num_atom))
    aa_residues   = []
    aa_chains     = []

    # ---------------
    # Local variables
    # ---------------
    i_resid        = 0
    curr_resid     = NaN
    curr_chain     = "?"
    curr_rname     = "    "
    residue_serial = NaN
    residue_name   = "    "
    chain_id       = "?"
    tmp_res_atoms  = []
    tmp_chain_res  = []
    segment_id     = " "

    # --------------------------------
    # Add atoms to residues and chains
    # --------------------------------
    for i_atom in 1:aa_num_atom + 1

        if i_atom <= aa_num_atom
            atom_name      = cif_data["_atom_site.label_atom_id"][i_atom]
            residue_name   = cif_data["_atom_site.label_comp_id"][i_atom]
            chain_id       = cif_data["_atom_site.label_asym_id"][i_atom]
            residue_serial = parse(Int, cif_data["_atom_site.label_seq_id"][i_atom])
            coor_x         = parse(Float64, cif_data["_atom_site.Cartn_x"][i_atom])
            coor_y         = parse(Float64, cif_data["_atom_site.Cartn_y"][i_atom])
            coor_z         = parse(Float64, cif_data["_atom_site.Cartn_z"][i_atom])
            segment_id     = " "
        end

        if chain_id != curr_chain || residue_serial > curr_resid + 1 || i_atom > aa_num_atom
            if length(tmp_res_atoms) > 0
                push!(aa_residues, AAResidue(curr_rname, tmp_res_atoms))
                tmp_res_atoms = []
            end
            if length(tmp_chain_res) > 0

                # -------------------------------
                # Determine chain molecule type
                # -------------------------------
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
            curr_chain = chain_id
            if i_atom > aa_num_atom
                break
            end
        end

        aa_atom_name[i_atom] = atom_name
        aa_coor[1, i_atom]   = coor_x
        aa_coor[2, i_atom]   = coor_y
        aa_coor[3, i_atom]   = coor_z

        if residue_serial != curr_resid
            i_resid += 1
            push!(tmp_chain_res, i_resid)
            if length(tmp_res_atoms) > 0
                push!(aa_residues, AAResidue(curr_rname, tmp_res_atoms))
                tmp_res_atoms = []
            end
            curr_resid = residue_serial
            curr_rname = residue_name
        end

        push!(tmp_res_atoms, i_atom)
    end

    new_molecule = AAMolecule(aa_atom_name, aa_coor, aa_residues, aa_chains)

    return new_molecule

end


# ===============
# Output CG mmCIF
# ===============
function write_mmCIF(top::GenTopology, conf::Conformation, system_name::AbstractString, args::Dict{String, <:Any}=Dict{String, Any}())

    verbose  = get(args, "verbose", false)

    cif_name = system_name * ".cif"
    cif_file = open(cif_name, "w")

    # ==========
    # write head
    # ==========
    println(cif_file, "data_$system_name")

    # ============
    # common lines
    # ============
    atom_site_dict_string = """#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
"""
    print(cif_file, atom_site_dict_string)


    # ====================
    # write CIF atom sites
    # ====================
    num_particles  = conf.num_particle

    chain_id_set  = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    tmp_chain_id  = 0
    tmp_seg_name  = ""
    real_chain_id = 1
    asym_id       = ""
    for i_bead in 1 : num_particles
        i_chain = top.top_atoms[i_bead].chain_id
        i_sname = top.top_atoms[i_bead].seg_name
        if i_chain != tmp_chain_id || i_sname != tmp_seg_name
            if tmp_chain_id > 0
                real_chain_id += 1
            end
            tmp_chain_id = i_chain
            tmp_seg_name = i_sname

            # ------------------
            # determine chain id
            # ------------------
            asym_id       = ""
            i_tmp = real_chain_id
            while i_tmp > 0
                (i_tmp, j) = divrem(i_tmp - 1, 52)
                asym_id *= chain_id_set[j + 1]
            end
        end

        @printf(cif_file,
                "ATOM %10d %2s %4s %1s %4s %4s %2d %8d %1s %15.3f %15.3f %15.3f %8.3f %8.3f %8.3f \n",
                i_bead,
                # top.top_atoms[i_bead].atom_type[1],
                "C",
                top.top_atoms[i_bead].atom_name,
                ".",
                top.top_atoms[i_bead].residue_name,
                asym_id,
                1,
                top.top_atoms[i_bead].residue_index,
                "?",
                conf.coors[1 , i_bead],
                conf.coors[2 , i_bead],
                conf.coors[3 , i_bead],
                0.0,
                0.0,
                top.top_atoms[i_bead].charge
                )
    end

    print(cif_file,"#\n")
    print(cif_file,"\n")

    close(cif_file)

    if verbose
        println(">           ... .pdb (CG) : DONE!")
    end

end
