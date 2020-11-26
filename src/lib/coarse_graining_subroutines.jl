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

function get_DNA3SPN_bond_length(bond_type::String, base_step::String)
    # Sugar-Base
    SB_length = Dict("A " => 4.864, "C " => 4.300, "G " => 4.973, "T " => 4.379)
    # Sugar-Phosphate
    SP_length = Dict(
        "AA" => 3.688, "AC" => 3.018, "AG" => 3.836, "AT" => 3.287,
        "CA" => 4.386, "CC" => 3.538, "CG" => 4.676, "CT" => 3.999,
        "GA" => 3.736, "GC" => 3.256, "GG" => 3.633, "GT" => 3.285,
        "TA" => 4.191, "TC" => 3.707, "TG" => 4.391, "TT" => 3.868
    )
    # Phosphate-Sugar
    PS_length = Dict(
        "XA" => 3.747, "XC" => 3.725, "XG" => 3.723, "XT" => 3.758,
        "AA" => 3.745, "AC" => 3.704, "AG" => 3.725, "AT" => 3.729,
        "CA" => 3.753, "CC" => 3.786, "CG" => 3.686, "CT" => 3.784,
        "GA" => 3.740, "GC" => 3.700, "GG" => 3.766, "GT" => 3.760,
        "TA" => 3.751, "TC" => 3.710, "TG" => 3.716, "TT" => 3.759
    )
    bond_length_data = Dict("SB"=>SB_length, "SP"=>SP_length, "PS"=>PS_length)

    return bond_length_data[bond_type][base_step]
end

function get_DNA3SPN_angle_equilibrium(angle_type::String, base_step::String)
    # Base-Sugar-Phosphate
    BSP_params = Dict(
        "AA" => 113.855, "AC" => 114.226, "AG" => 112.201, "AT" => 111.931,
        "CA" => 113.822, "CC" => 112.056, "CG" => 116.081, "CT" => 111.008,
        "GA" => 114.665, "GC" => 118.269, "GG" => 110.102, "GT" => 111.146,
        "TA" => 113.984, "TC" => 115.457, "TG" => 113.397, "TT" => 113.606
    )
    # Phosphate-Sugar-Base
    PSB_params = Dict(
        "XA" => 108.200, "XC" => 103.850, "XG" => 111.750, "XT" => 98.523,
        "AA" => 108.826, "AC" => 105.066, "AG" => 112.796, "AT" => 99.442,
        "CA" => 107.531, "CC" => 103.509, "CG" => 110.594, "CT" => 97.807,
        "GA" => 108.064, "GC" => 103.135, "GG" => 112.654, "GT" => 98.577,
        "TA" => 108.414, "TC" => 103.853, "TG" => 111.732, "TT" => 98.271
    )
    # Phosphate-Sugar-Phosphate
    PSP_params = Dict(
        # TODO: currently using "X" as "A", should be changed to average???
        "XAA" => 120.685, "XAC" => 112.882, "XAG" => 113.827, "XAT" => 117.435,
        "XCA" => 119.061, "XCC" => 120.353, "XCG" => 113.240, "XCT" => 121.103,
        "XGA" => 122.182, "XGC" => 118.658, "XGG" => 120.489, "XGT" => 122.928,
        "XTA" => 117.235, "XTC" => 112.084, "XTG" => 111.714, "XTT" => 119.324,
        "AAA" => 120.685, "AAC" => 112.882, "AAG" => 113.827, "AAT" => 117.435,
        "ACA" => 119.061, "ACC" => 120.353, "ACG" => 113.240, "ACT" => 121.103,
        "AGA" => 122.182, "AGC" => 118.658, "AGG" => 120.489, "AGT" => 122.928,
        "ATA" => 117.235, "ATC" => 112.084, "ATG" => 111.714, "ATT" => 119.324,
        "CAA" => 122.866, "CAC" => 115.083, "CAG" => 116.036, "CAT" => 119.640,
        "CCA" => 120.442, "CCC" => 121.712, "CCG" => 114.602, "CCT" => 122.446,
        "CGA" => 124.721, "CGC" => 121.204, "CGG" => 122.937, "CGT" => 125.429,
        "CTA" => 119.317, "CTC" => 114.156, "CTG" => 113.756, "CTT" => 121.413,
        "GAA" => 120.809, "GAC" => 112.897, "GAG" => 113.816, "GAT" => 117.461,
        "GCA" => 119.550, "GCC" => 120.788, "GCG" => 113.687, "GCT" => 121.506,
        "GGA" => 121.512, "GGC" => 118.019, "GGG" => 119.634, "GGT" => 122.157,
        "GTA" => 117.087, "GTC" => 111.922, "GTG" => 111.501, "GTT" => 119.185,
        "TAA" => 122.361, "TAC" => 114.671, "TAG" => 115.653, "TAT" => 119.219,
        "TCA" => 121.235, "TCC" => 122.532, "TCG" => 115.417, "TCT" => 123.284,
        "TGA" => 123.936, "TGC" => 120.395, "TGG" => 122.319, "TGT" => 124.730,
        "TTA" => 119.004, "TTC" => 113.847, "TTG" => 113.465, "TTT" => 121.093
    )
    # Sugar-Phosphate-Sugar
    SPS_params = Dict(
        "AA" => 94.805, "AC" => 94.462, "AG" => 95.308, "AT" => 95.232,
        "CA" => 95.110, "CC" => 98.906, "CG" => 92.244, "CT" => 97.476,
        "GA" => 94.973, "GC" => 92.666, "GG" => 97.929, "GT" => 97.640,
        "TA" => 94.886, "TC" => 93.066, "TG" => 93.999, "TT" => 95.122
    )
    angle_data = Dict(
        "BSP" => BSP_params,
        "PSB" => PSB_params,
        "PSP" => PSP_params,
        "SPS" => SPS_params
    )

    return angle_data[angle_type][base_step]
end

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
        "XA" => 292, "XT" => 407, "XC" => 359, "XG" => 280,
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

function get_DNA3SPN_dihedral_equilibrium(angle_type::String, base_step::String)
    # Base-Sugar-Phosphate-Sugar
    BSPS_params = Dict(
        "AA" => -23.347, "AC" => -27.858, "AG" => -27.117, "AT" => -29.246,
        "CA" => -31.608, "CC" => -31.364, "CG" => -34.383, "CT" => -33.819,
        "GA" => -16.641, "GC" => -17.077, "GG" => -20.529, "GT" => -21.472,
        "TA" => -36.960, "TC" => -39.034, "TG" => -39.283, "TT" => -38.799
    )
    # Sugar-Phosphate-Sugar-Base
    SPSB_params = Dict(
        "AA" => 45.425, "AC" => 54.789, "AG" => 46.984, "AT" => 57.208,
        "CA" => 45.195, "CC" => 49.771, "CG" => 44.547, "CT" => 53.367,
        "GA" => 41.089, "GC" => 45.515, "GG" => 43.923, "GT" => 51.560,
        "TA" => 47.078, "TC" => 52.838, "TG" => 46.053, "TT" => 54.408
    )
    # Sugar-Phosphate-Sugar-Phosphate
    SPSP_params = Dict(
        "AAA" =>  179.785, "AAC" =>  173.331, "AAG" =>  171.377, "AAT" =>  173.860,
        "ACA" => -176.300, "ACC" => -177.745, "ACG" => -177.543, "ACT" => -178.626,
        "AGA" => -169.949, "AGC" => -168.414, "AGG" =>  179.834, "AGT" => -175.422,
        "ATA" => -179.491, "ATC" =>  179.733, "ATG" =>  177.177, "ATT" => -178.801,
        "CAA" => -179.648, "CAC" =>  173.730, "CAG" =>  171.730, "CAT" =>  174.273,
        "CCA" =>  178.306, "CCC" =>  176.814, "CCG" =>  177.164, "CCT" =>  175.901,
        "CGA" => -172.058, "CGC" => -170.458, "CGG" =>  177.459, "CGT" => -177.699,
        "CTA" =>  176.566, "CTC" =>  175.846, "CTG" =>  173.254, "CTT" =>  177.252,
        "GAA" =>  174.706, "GAC" =>  168.411, "GAG" =>  166.426, "GAT" =>  168.841,
        "GCA" =>  173.035, "GCC" =>  171.523, "GCG" =>  172.053, "GCT" =>  170.601,
        "GGA" => -174.234, "GGC" => -172.619, "GGG" =>  175.674, "GGT" => -179.682,
        "GTA" =>  174.167, "GTC" =>  173.514, "GTG" =>  170.969, "GTT" =>  174.788,
        "TAA" => -177.232, "TAC" =>  176.044, "TAG" =>  174.067, "TAT" =>  176.665,
        "TCA" => -177.663, "TCC" => -179.135, "TCG" => -178.953, "TCT" =>  179.965,
        "TGA" => -169.881, "TGC" => -168.369, "TGG" =>  179.680, "TGT" => -175.458,
        "TTA" =>  177.790, "TTC" =>  177.037, "TTG" =>  174.444, "TTT" =>  178.479
    )
    # Phosphate-Sugar-Phosphate-Sugar
    PSPS_params = Dict(
        # TODO: currently using "X" as "A", should be changed to average???
        "XAA" => -155.622, "XAC" => -152.885, "XAG" => -151.259, "XAT" => -156.185,
        "XCA" => -156.388, "XCC" => -155.577, "XCG" => -156.063, "XCT" => -157.660,
        "XGA" => -159.083, "XGC" => -159.751, "XGG" => -154.497, "XGT" => -159.668,
        "XTA" => -152.487, "XTC" => -151.938, "XTG" => -150.672, "XTT" => -155.597,
        "AAA" => -155.622, "AAC" => -152.885, "AAG" => -151.259, "AAT" => -156.185,
        "ACA" => -156.388, "ACC" => -155.577, "ACG" => -156.063, "ACT" => -157.660,
        "AGA" => -159.083, "AGC" => -159.751, "AGG" => -154.497, "AGT" => -159.668,
        "ATA" => -152.487, "ATC" => -151.938, "ATG" => -150.672, "ATT" => -155.597,
        "CAA" => -156.021, "CAC" => -152.981, "CAG" => -151.273, "CAT" => -156.309,
        "CCA" => -155.364, "CCC" => -154.499, "CCG" => -155.058, "CCT" => -156.547,
        "CGA" => -158.746, "CGC" => -159.509, "CGG" => -153.638, "CGT" => -159.033,
        "CTA" => -151.817, "CTC" => -151.269, "CTG" => -149.902, "CTT" => -154.955,
        "GAA" => -154.534, "GAC" => -151.854, "GAG" => -150.223, "GAT" => -155.116,
        "GCA" => -154.009, "GCC" => -153.155, "GCG" => -153.791, "GCT" => -155.211,
        "GGA" => -157.783, "GGC" => -158.478, "GGG" => -153.379, "GGT" => -158.439,
        "GTA" => -151.220, "GTC" => -150.726, "GTG" => -149.471, "GTT" => -154.288,
        "TAA" => -156.903, "TAC" => -153.864, "TAG" => -152.178, "TAT" => -157.225,
        "TCA" => -156.627, "TCC" => -155.754, "TCG" => -156.236, "TCT" => -157.799,
        "TGA" => -159.780, "TGC" => -160.478, "TGG" => -154.803, "TGT" => -160.164,
        "TTA" => -152.217, "TTC" => -151.655, "TTG" => -150.303, "TTT" => -155.342
    )
    angle_data = Dict(
        "BSPS" => BSPS_params,
        "SPSB" => SPSB_params,
        "SPSP" => SPSP_params,
        "PSPS" => PSPS_params
    )

    return angle_data[angle_type][base_step]

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

# ==============
# PWMcos contact
# ==============

function is_protein_DNA_Go_contact(resid1_atoms::Vector{Int}, resid2_atoms::Vector{Int}, atom_names::Vector{String}, atom_coors::Array{<:Real, 2})
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
            if dist_12 < pro_DNA_GO_ATOMIC_CUTOFF
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




