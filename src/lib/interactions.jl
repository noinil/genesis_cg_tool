###############################################################################
#                        Functions related to CG models                       #
###############################################################################


###############################################################################
#                       ____            _       _                             #
#                      |  _ \ _ __ ___ | |_ ___(_)_ __                        #
#                      | |_) | '__/ _ \| __/ _ \ | '_ \                       #
#                      |  __/| | | (_) | ||  __/ | | | |                      #
#                      |_|   |_|  \___/ \__\___|_|_| |_|                      #
#                                                                             #
###############################################################################

# ------------------------------
# General Protein Native Contact
# ------------------------------

function is_protein_native_contact(resid1_atoms::Vector{Int}, resid2_atoms::Vector{Int}, atom_names::Vector{String}, atom_coors::Array{<:Real, 2})
    for i in resid1_atoms
        atom_name_1 = atom_names[i]
        if atom_name_1[1] == 'H'
            continue
        end
        coor_1 = atom_coors[:, i]
        for j in resid2_atoms
            atom_name_2 = atom_names[j]
            if atom_name_2[1] == 'H'
                continue
            end
            coor_2  = atom_coors[:, j]
            dist_12 = compute_distance(coor_1, coor_2)
            if dist_12 < AICG_GO_ATOMIC_CUTOFF
                return true
            end
        end
    end
    return false
end


# --------------------
# AICG2+ Protein Model
# --------------------

function is_protein_backbone(atom_name::String)
    if in(atom_name, ("N", "C", "O", "OXT", "CA"))
        return true
    end
    return false
end

function is_protein_hb_donor(atom_name::String, res_name::String)
    if atom_name[1] == 'N'
        return true
    elseif atom_name[1] == 'S' && res_name == "CYS"
        return true
    elseif atom_name[1] == 'O'
        if  ( res_name == "SER" && atom_name == "OG"  ) ||
            ( res_name == "THR" && atom_name == "OG1" ) ||
            ( res_name == "TYR" && atom_name == "OH"  )
            return true
        end
    end
    return false
end

function is_protein_hb_acceptor(atom_name::String)
    if atom_name[1] == 'O' || atom_name[1] == 'S'
        return true
    end
    return false
end

function is_protein_cation(atom_name::String, res_name::String)
    if atom_name[1] == 'N'
        if  ( res_name == "ARG" && atom_name == "NH1" ) ||
            ( res_name == "ARG" && atom_name == "NH2" ) ||
            ( res_name == "LYS" && atom_name == "NZ"  )
            return true
        end
    end
    return false
end

function is_protein_anion(atom_name::String, res_name::String)
    if atom_name[1] == 'O'
        if  ( res_name == "GLU" && atom_name == "OE1" ) ||
            ( res_name == "GLU" && atom_name == "OE2" ) ||
            ( res_name == "ASP" && atom_name == "OD1" ) ||
            ( res_name == "ASP" && atom_name == "OD2" )
            return true
        end
    end
    return false
end

function is_protein_hb_pair(atom_name_1::String, res_name_1::String, atom_name_2::String, res_name_2::String)
    if  is_protein_hb_acceptor(atom_name_1) &&
        is_protein_hb_donor(atom_name_2, res_name_2)
        return true
    elseif is_protein_hb_acceptor(atom_name_2) &&
        is_protein_hb_donor(atom_name_1, res_name_1)
        return true
    end
    return false
end

function is_protein_sb_pair(atom_name_1::String, res_name_1::String, atom_name_2::String, res_name_2::String)
    if  is_protein_cation(atom_name_1, res_name_1) &&
        is_protein_anion(atom_name_2,  res_name_2)
        return true
    elseif is_protein_cation(atom_name_2, res_name_2) &&
        is_protein_anion(atom_name_1,  res_name_1)
        return true
    end
    return false
end

function is_protein_nonsb_charge_pair(atom_name_1::String, res_name_1::String, atom_name_2::String, res_name_2::String)
    if  is_protein_cation(atom_name_1, res_name_1) ||
        is_protein_anion(atom_name_1,  res_name_1) ||
        is_protein_cation(atom_name_2, res_name_2) ||
        is_protein_anion(atom_name_2,  res_name_2)
        return true
    end
    return false
end

function count_aicg_atomic_contact(resid1_atoms::Vector{Int}, resid2_atoms::Vector{Int}, res_name_1::String, res_name_2::String, atom_names::Vector{String}, atom_coors::Array{<:Real, 2})
    contact_count                   = zeros(Int, 17)
    contact_count[AICG_ITYPE_OFFST] = 1
    num_short_range_contact         = 0
    for i in resid1_atoms
        atom_name_1 = atom_names[i]
        if atom_name_1[1] == 'H'
            continue
        end
        coor_1 = atom_coors[:, i]
        for j in resid2_atoms
            atom_name_2     = atom_names[j]
            if atom_name_2[1] == 'H'
                continue
            end
            coor_2          = atom_coors[:, j]
            dist_12         = compute_distance(coor_1, coor_2)

            is_hb           = is_protein_hb_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
            is_sb           = is_protein_sb_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
            is_nonsb_charge = is_protein_nonsb_charge_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
            is_1_backbone   = is_protein_backbone(atom_name_1)
            is_2_backbone   = is_protein_backbone(atom_name_2)
            if dist_12 < AICG_GO_ATOMIC_CUTOFF
                contact_count[AICG_ITYPE_LR_CT] += 1
            end
            if dist_12 < AICG_ATOMIC_CUTOFF
                num_short_range_contact += 1
                if is_1_backbone && is_2_backbone
                    if is_hb
                        if dist_12 < AICG_HYDROGEN_BOND_CUTOFF
                            contact_count[AICG_ITYPE_BB_HB] += 1
                        else
                            contact_count[AICG_ITYPE_BB_DA] += 1
                        end
                    elseif atom_name_1[1] == 'C' || atom_name_2[1] == 'C'
                        contact_count[AICG_ITYPE_BB_CX] += 1
                    else
                        contact_count[AICG_ITYPE_BB_XX] += 1
                    end
                elseif ( !is_1_backbone ) && ( !is_2_backbone )
                    if is_hb
                        if is_sb
                            if dist_12 < AICG_SALT_BRIDGE_CUTOFF
                                contact_count[AICG_ITYPE_SS_SB] += 1
                            else
                                contact_count[AICG_ITYPE_SS_QX] += 1
                            end
                        elseif dist_12 < AICG_HYDROGEN_BOND_CUTOFF
                            contact_count[AICG_ITYPE_SS_HB] += 1
                        elseif is_nonsb_charge
                            contact_count[AICG_ITYPE_SS_QX] += 1
                        else
                            contact_count[AICG_ITYPE_SS_DA] += 1
                        end
                    elseif is_nonsb_charge
                        contact_count[AICG_ITYPE_SS_QX] += 1
                    elseif atom_name_1[1] == 'C' || atom_name_2[1] == 'C'
                        contact_count[AICG_ITYPE_SS_CX] += 1
                    else
                        contact_count[AICG_ITYPE_SS_XX] += 1
                    end
                elseif ( is_1_backbone && ( !is_2_backbone ) ) ||
                    ( is_2_backbone && ( !is_1_backbone ) )
                    if is_hb
                        if dist_12 < AICG_HYDROGEN_BOND_CUTOFF
                            contact_count[AICG_ITYPE_SB_HB] += 1
                        elseif is_nonsb_charge
                            contact_count[AICG_ITYPE_SB_QX] += 1
                        else
                            contact_count[AICG_ITYPE_SB_DA] += 1
                        end
                    elseif is_nonsb_charge
                        contact_count[AICG_ITYPE_SB_QX] += 1
                    elseif atom_name_1[1] == 'C' || atom_name_2[1] == 'C'
                        contact_count[AICG_ITYPE_SB_CX] += 1
                    else
                        contact_count[AICG_ITYPE_SB_XX] += 1
                    end
                end
            end
        end
    end

    # control the number of long-range contacts
    if AICG_GO_ATOMIC_CUTOFF > AICG_ATOMIC_CUTOFF
        contact_count[AICG_ITYPE_LR_CT] -= num_short_range_contact
    else
        contact_count[AICG_ITYPE_LR_CT]  = 0
    end

    # control the number of salty bridge
    if contact_count[AICG_ITYPE_SS_SB]  >= 2
        contact_count[AICG_ITYPE_SS_QX] += contact_count[AICG_ITYPE_SS_SB] - 1
        contact_count[AICG_ITYPE_SS_SB]  = 1
    end

    return contact_count
end


###############################################################################
#                           ____   _   _     _                                #
#                          |  _ \ | \ | |   / \                               #
#                          | | | ||  \| |  / _ \                              #
#                          | |_| || |\  | / ___ \                             #
#                          |____/ |_| \_|/_/   \_\                            #
#                                                                             #
###############################################################################

# -----------------
# 3SPN.2C DNA model
# -----------------

function get_DNA3SPN_angle_param(angle_type::String, base_step::String)
    # Base-Sugar-Phosphate
    BSP_params = Dict(
        "AA" => 460, "AT" => 370, "AC" => 442, "AG" => 358,
        "TA" => 120, "TT" => 460, "TC" => 383, "TG" => 206,
        "CA" => 206, "CT" => 358, "CC" => 278, "CG" => 278,
        "GA" => 383, "GT" => 442, "GC" => 336, "GG" => 278
    )
    # Phosphate-Sugar-Base
    PSB_params = Dict(
        "AA" => 460, "TA" => 120, "CA" => 206, "GA" => 383,
        "AT" => 370, "TT" => 460, "CT" => 358, "GT" => 442,
        "AC" => 442, "TC" => 383, "CC" => 278, "GC" => 336,
        "AG" => 358, "TG" => 206, "CG" => 278, "GG" => 278
    )
    # Phosphate-Sugar-Phosphate
    PSP_params = Dict(
        "all" => 300
    )
    # Sugar-Phosphate-Sugar
    SPS_params = Dict(
        "AA" => 355, "AT" => 147, "AC" => 464, "AG" => 368,
        "TA" => 230, "TT" => 355, "TC" => 442, "TG" => 273,
        "CA" => 273, "CT" => 368, "CC" => 165, "CG" => 478,
        "GA" => 442, "GT" => 464, "GC" => 228, "GG" => 165
    )
    angle_params = Dict(
        "BSP" => BSP_params,
        "PSB" => PSB_params,
        "PSP" => PSP_params,
        "SPS" => SPS_params
    )

    return angle_params[angle_type][base_step] * JOU2CAL
end


###############################################################################
#                              ____   _   _     _                             #
#                             |  _ \ | \ | |   / \                            #
#                             | |_) ||  \| |  / _ \                           #
#                             |  _ < | |\  | / ___ \                          #
#                             |_| \_\|_| \_|/_/   \_\                         #
#                                                                             #
###############################################################################

# -------------------------
# RNA structure-based model
# -------------------------
function is_RNA_hydrogen_bond(atom_name_1::Char, atom_name_2::Char)
    special_atom_list = ['F', 'O', 'N']
    if atom_name_1 in special_atom_list && atom_name_2 in special_atom_list
        return true
    end
    return false
end

function compute_RNA_native_contact(resid1_atoms::Vector{Int}, resid2_atoms::Vector{Int}, atom_names::Vector{String}, atom_coors::Array{<:Real, 2})
    hb_count = 0
    min_dist = 1.0e50
    for i in resid1_atoms
        atom_name_1 = atom_names[i]
        if atom_name_1[1] == 'H'
            continue
        end
        coor_1 = atom_coors[:, i]
        for j in resid2_atoms
            atom_name_2 = atom_names[j]
            if atom_name_2[1] == 'H'
                continue
            end
            coor_2 = atom_coors[:, j]
            dist_12 = compute_distance(coor_1, coor_2)
            if dist_12 < RNA_GO_ATOMIC_CUTOFF && is_RNA_hydrogen_bond(atom_name_1[1], atom_name_2[1])
                hb_count += 1
            end
            if dist_12 < min_dist
                min_dist = dist_12
            end
        end
    end
    return (min_dist, hb_count)
end


###############################################################################
#                        _         _                ____   _   _     _        #
#     _ __   _ __  ___  | |_  ___ (_) _ __         |  _ \ | \ | |   / \       #
#    | '_ \ | '__|/ _ \ | __|/ _ \| || '_ \  _____ | | | ||  \| |  / _ \      #
#    | |_) || |  | (_) || |_|  __/| || | | ||_____|| |_| || |\  | / ___ \     #
#    | .__/ |_|   \___/  \__|\___||_||_| |_|       |____/ |_| \_|/_/   \_\    #
#    |_|                                                                      #
###############################################################################

# ==============
# PWMcos contact
# ==============

function is_PWMcos_contact(resid1_atoms::Vector{Int}, resid2_atoms::Vector{Int}, atom_names::Vector{String}, atom_coors::Array{<:Real, 2})
    for i in resid1_atoms
        atom_name_1 = atom_names[i]
        if atom_name_1[1] == 'H'
            continue
        end
        coor_1 = atom_coors[:, i]
        for j in resid2_atoms
            atom_name_2 = atom_names[j]
            if atom_name_2[1] == 'H'
                continue
            end
            coor_2  = atom_coors[:, j]
            dist_12 = compute_distance(coor_1, coor_2)
            if dist_12 < PWMCOS_ATOMIC_CUTOFF
                return true
            end
        end
    end
    return false
end


###############################################################################
#                        _         _                ____   _   _     _        #
#     _ __   _ __  ___  | |_  ___ (_) _ __         |  _ \ | \ | |   / \       #
#    | '_ \ | '__|/ _ \ | __|/ _ \| || '_ \  _____ | |_) ||  \| |  / _ \      #
#    | |_) || |  | (_) || |_|  __/| || | | ||_____||  _ < | |\  | / ___ \     #
#    | .__/ |_|   \___/  \__|\___||_||_| |_|       |_| \_\|_| \_|/_/   \_\    #
#    |_|                                                                      #
###############################################################################

function is_protein_RNA_native_contact(resid1_atoms::Vector{Int}, resid2_atoms::Vector{Int}, atom_names::Vector{String}, atom_coors::Array{<:Real, 2})
    for i in resid1_atoms
        atom_name_1 = atom_names[i]
        if atom_name_1[1] == 'H'
            continue
        end
        coor_1 = atom_coors[:, i]
        for j in resid2_atoms
            atom_name_2 = atom_names[j]
            if atom_name_2[1] == 'H'
                continue
            end
            coor_2  = atom_coors[:, j]
            dist_12 = compute_distance(coor_1, coor_2)
            if dist_12 < AICG_GO_ATOMIC_CUTOFF
                return true
            end
        end
    end
    return false
end




