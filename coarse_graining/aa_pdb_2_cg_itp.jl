#!/usr/bin/env julia

using ArgParse
using LinearAlgebra

###########################################################################
#                          Force Field Parameters                         #
###########################################################################

# ==================
# Physical Constants
# ==================
const CAL2JOU = 4.184

# =====================================
# General Parameters: Mass, Charge, ...
# =====================================

res_mass_dict = Dict(
    "ALA" =>  71.09,
    "ARG" => 156.19,
    "ASN" => 114.11,
    "ASP" => 115.09,
    "CYS" => 103.15,
    "GLN" => 128.14,
    "GLU" => 129.12,
    "GLY" =>  57.05,
    "HIS" => 137.14,
    "ILE" => 113.16,
    "LEU" => 113.16,
    "LYS" => 128.17,
    "MET" => 131.19,
    "PHE" => 147.18,
    "PRO" =>  97.12,
    "SER" =>  87.08,
    "THR" => 101.11,
    "TRP" => 186.21,
    "TYR" => 163.18,
    "VAL" =>  99.14
)

res_charge_dict = Dict(
    "ALA" =>  0.0,
    "ARG" =>  1.0,
    "ASN" =>  0.0,
    "ASP" => -1.0,
    "CYS" =>  0.0,
    "GLN" =>  0.0,
    "GLU" => -1.0,
    "GLY" =>  0.0,
    "HIS" =>  0.0,
    "ILE" =>  0.0,
    "LEU" =>  0.0,
    "LYS" =>  1.0,
    "MET" =>  0.0,
    "PHE" =>  0.0,
    "PRO" =>  0.0,
    "SER" =>  0.0,
    "THR" =>  0.0,
    "TRP" =>  0.0,
    "TYR" =>  0.0,
    "VAL" =>  0.0
)

# ===============================
# Protein AICG2+ Model Parameters
# ===============================

# AICG2+ bond force constant
const aicg_bond_k               = 110.40 * CAL2JOU * 100 * 2
# AICG2+ sigma for Gaussian angle
const aicg_ang_gauss_sigma      = 0.15 * 0.1  # nm
# AICG2+ sigma for Gaussian dihedral
const aicg_dih_gauss_sigma      = 0.15        # Rad ??
# AICG2+ atomistic contact cutoff
const aicg_go_atomic_cutoff     = 6.5
# AICG2+ pairwise interaction cutoff
const aicg_atomic_cutoff        = 5.0
# AICG2+ hydrogen bond cutoff
const aicg_hydrogen_bond_cutoff = 3.2
# AICG2+ salt bridge cutoff
const aicg_salt_bridge_cutoff   = 3.5
# AICG2+ energy cutoffs
const aicg_ene_upper_lim        = -0.5
const aicg_ene_lower_lim        = -5.0
# average and general AICG2+ energy values
const aicg_13_ave               = 1.72
const aicg_14_ave               = 1.23
const aicg_contact_ave          = 0.55
const aicg_13_gen               = 1.11
const aicg_14_gen               = 0.87
const aicg_contact_gen          = 0.32

# AICG2+ pairwise interaction pairs
const aicg_itype_bb_hb = 1  # B-B hydrogen bonds
const aicg_itype_bb_da = 2  # B-B donor-accetor contacts
const aicg_itype_bb_cx = 3  # B-B carbon-X contacts
const aicg_itype_bb_xx = 4  # B-B other
const aicg_itype_ss_hb = 5  # S-S hydrogen bonds
const aicg_itype_ss_sb = 6  # S-S salty bridge
const aicg_itype_ss_da = 7  # S-S donor-accetor contacts
const aicg_itype_ss_cx = 8  # S-S carbon-X contacts
const aicg_itype_ss_qx = 9  # S-S charge-X contacts
const aicg_itype_ss_xx = 10 # S-S other
const aicg_itype_sb_hb = 11 # S-B hydrogen bonds
const aicg_itype_sb_da = 12 # S-B donor-accetor contacts
const aicg_itype_sb_cx = 13 # S-B carbon-X contacts
const aicg_itype_sb_qx = 14 # S-B charge-X contacts
const aicg_itype_sb_xx = 15 # S-B other
const aicg_itype_lr_ct = 16 # long range contacts
const aicg_itype_offst = 17  # offset

aicg_pairwise_energy = ones(Float64, (17, 1))
aicg_pairwise_energy[aicg_itype_bb_hb] = - 1.4247   # B-B hydrogen bonds
aicg_pairwise_energy[aicg_itype_bb_da] = - 0.4921   # B-B donor-accetor contacts
aicg_pairwise_energy[aicg_itype_bb_cx] = - 0.2404   # B-B carbon-X contacts
aicg_pairwise_energy[aicg_itype_bb_xx] = - 0.1035   # B-B other
aicg_pairwise_energy[aicg_itype_ss_hb] = - 5.7267   # S-S hydrogen bonds
aicg_pairwise_energy[aicg_itype_ss_sb] = -12.4878   # S-S salty bridge
aicg_pairwise_energy[aicg_itype_ss_da] = - 0.0308   # S-S donor-accetor contacts
aicg_pairwise_energy[aicg_itype_ss_cx] = - 0.1113   # S-S carbon-X contacts
aicg_pairwise_energy[aicg_itype_ss_qx] = - 0.2168   # S-S charge-X contacts
aicg_pairwise_energy[aicg_itype_ss_xx] =   0.2306   # S-S other
aicg_pairwise_energy[aicg_itype_sb_hb] = - 3.4819   # S-B hydrogen bonds
aicg_pairwise_energy[aicg_itype_sb_da] = - 0.1809   # S-B donor-accetor contacts
aicg_pairwise_energy[aicg_itype_sb_cx] = - 0.1209   # S-B carbon-X contacts
aicg_pairwise_energy[aicg_itype_sb_qx] = - 0.2984   # S-B charge-X contacts
aicg_pairwise_energy[aicg_itype_sb_xx] = - 0.0487   # S-B other
aicg_pairwise_energy[aicg_itype_lr_ct] = - 0.0395   # long range contacts
aicg_pairwise_energy[aicg_itype_offst] = - 0.1051   # offset


# ====================
# GRO TOP File Options
# ====================

# "NREXCL" in "[moleculetype]"
const mol_nr_excl            = 3
# "CGNR" in "[atoms]"
const aicg_atom_func_nr      = 1
# "f" in "[bonds]"
const aicg_bond_func_type    = 1
# "f" in AICG-type "[angles]"
const aicg_ang_g_func_type   = 21
# "f" in Flexible-type "[angles]"
const aicg_ang_f_func_type   = 22
# "f" in AICG-type "[dihedral]"
const aicg_dih_g_func_type   = 21
# "f" in Flexible-type "[dihedral]"
const aicg_dih_f_func_type   = 22
# "f" in Go-contacts "[pairs]"
const aicg_contact_func_type = 2


###############################################################################
#                                  Functions                                  #
###############################################################################

# =============================
# Parsing Commandline Arguments
# =============================
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "pdb"
            help     = "PDB file name."
            required = true
            arg_type = String
    end
    return parse_args(s)
end

# ===================
# Geometric Functions
# ===================

# --------
# Distance
# --------
function compute_distance(coor1, coor2)
    d = coor1 - coor2
    return norm(d)
end

# -----
# Angle
# -----
function compute_angle(coor1, coor2, coor3)
    v1 = coor1 - coor2
    v2 = coor3 - coor2
    n1 = norm(v1)
    n2 = norm(v2)
    return acos( dot(v1, v2) / n1 / n2) / pi * 180.0
end

# --------
# Dihedral
# --------
function compute_dihedral(coor1, coor2, coor3, coor4)
    v12   = coor2 - coor1
    v23   = coor3 - coor2
    v34   = coor4 - coor3
    c123  = cross(v12, v23)
    c234  = cross(v23, v34)
    nc123 = norm(c123)
    nc234 = norm(c234)
    dih   = acos( dot(c123, c234) / nc123 / nc234)
    c1234 = cross(c123, c234)
    judge = dot(c1234, v23)
    dih   = judge < 0 ? - dih : dih
    return dih / pi * 180.0
end

# ===============================
# Structural Biological Functions
# ===============================

# --------------------
# AICG2+ Protein Model
# --------------------
function is_protein_backbone(atom_name)
    if in(atom_name, ("N", "C", "O", "OXT", "CA"))
        return true
    end
    return false
end
function is_protein_hb_donor(atom_name, res_name)
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
    else
        return false
    end
end
function is_protein_hb_acceptor(atom_name)
    if atom_name[1] == 'O' || atom_name[1] == 'S'
        return true
    end
    return false
end
function is_protein_cation(atom_name, res_name)
    if atom_name[1] == 'N'
        if  ( res_name == "ARG" && atom_name == "NH1" ) || 
            ( res_name == "ARG" && atom_name == "NH2" ) || 
            ( res_name == "LYS" && atom_name == "NZ"  )
            return true
        end
    else
        return false
    end
end
function is_protein_anion(atom_name, res_name)
    if atom_name[1] == 'O'
        if  ( res_name == "GLU" && atom_name == "OE1" ) ||
            ( res_name == "GLU" && atom_name == "OE2" ) ||
            ( res_name == "ASP" && atom_name == "OD1" ) ||
            ( res_name == "ASP" && atom_name == "OD2" )
            return true
        end
    else
        return false
    end
end
function is_protein_hb_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
    if  is_protein_hb_acceptor(atom_name_1) &&
        is_protein_hb_donor(atom_name_2, res_name_2)
        return true
    elseif is_protein_hb_acceptor(atom_name_2) &&
        is_protein_hb_donor(atom_name_1, res_name_1)
        return true
    else
        return false
    end
end
function is_protein_sb_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
    if  is_protein_cation(atom_name_1, res_name_1) &&
        is_protein_anion(atom_name_2,  res_name_2)
        return true
    elseif is_protein_cation(atom_name_2, res_name_2) &&
        is_protein_anion(atom_name_1,  res_name_1)
        return true
    else
        return false
    end
end
function is_protein_nonsb_charge_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
    if  is_protein_cation(atom_name_1, res_name_1) ||
        is_protein_anion(atom_name_1,  res_name_1) ||
        is_protein_cation(atom_name_2, res_name_2) ||
        is_protein_anion(atom_name_2,  res_name_2)
        return True
    else
        return False
    end
end
function is_protein_go_contact(resid1, resid2, atom_names, atom_coors)
    for i in resid1.first : resid1.last
        atom_name_1 = atom_names[i]
        if atom_name_1[1] == 'H'
            continue
        end
        coor_1 = atom_coors[i]
        for j in resid2.fist : resid2.last
            atom_name_2 = atom_names[j]
            if atom_name_2[1] == 'H'
                continue
            end
            coor_2  = atom_coors[j]
            dist_12 = compute_distance(coor_1, coor_2)
            if dist_12 < aicg_go_atomic_cutoff
                return True
            end
        end
    end
    return False
end
function count_aicg_atomic_contact(resid1, resid2, res_name_1, res_name_2, atom_names, atom_coors)
    contact_count                   = [0 for _ in 1:17]
    contact_count[aicg_itype_offst] = 1
    num_short_range_contact         = 0
    for i in resid1.first : resid1.last
        atom_name_1 = atom_names[i]
        if atom_name_1[1] == 'H'
            continue
        end
        coor_1 = atom_coors[i]
        for j in resid2.start : resid2.last
            atom_name_2     = atom_names[j]
            if atom_name_2[1] == 'H'
                continue
            end
            coor_2          = atom_coors[j]
            dist_12         = compute_distance(coor_1, coor_2)

            is_hb           = is_protein_hb_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
            is_sb           = is_protein_sb_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
            is_nonsb_charge = is_protein_nonsb_charge_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
            is_1_backbone   = is_protein_backbone(atom_name_1)
            is_2_backbone   = is_protein_backbone(atom_name_2)
            if dist_12 < aicg_go_atomic_cutoff
                contact_count[aicg_itype_lr_ct] += 1
            elseif dist_12 < aicg_atomic_cutoff
                num_short_range_contact += 1
                if is_1_backbone && is_2_backbone
                    if is_hb
                        if dist_12 < aicg_hydrogen_bond_cutoff
                            contact_count[aicg_itype_bb_hb] += 1
                        else
                            contact_count[aicg_itype_bb_da] += 1
                        end
                    elseif atom_name_1[1] == 'C' || atom_name_2[1] == 'C'
                        contact_count[aicg_itype_bb_cx] += 1
                    else
                        contact_count[aicg_itype_bb_xx] += 1
                    end
                elseif ( !is_1_backbone ) && ( !is_2_backbone )
                    if is_hb
                        if is_sb
                            if dist_12 < aicg_salt_bridge_cutoff
                                contact_count[aicg_itype_ss_sb] += 1
                            else
                                contact_count[aicg_itype_ss_qx] += 1
                            end
                        elseif dist_12 < aicg_hydrogen_bond_cutoff
                            contact_count[aicg_itype_ss_hb] += 1
                        elseif is_nonsb_charge
                            contact_count[aicg_itype_ss_qx] += 1
                        else
                            contact_count[aicg_itype_ss_da] += 1
                        end
                    elseif is_nonsb_charge
                        contact_count[aicg_itype_ss_qx] += 1
                    elseif atom_name_1[1] == 'C' || atom_name_2[1] == 'C'
                        contact_count[aicg_itype_ss_cx] += 1
                    else
                        contact_count[aicg_itype_ss_xx] += 1
                    end
                elseif ( is_1_backbone && ( !is_2_backbone ) ) ||
                    ( is_2_backbone && ( !is_1_backbone ) )
                    if is_hb
                        if dist_12 < aicg_hydrogen_bond_cutoff
                            contact_count[aicg_itype_sb_hb] += 1
                        elseif is_nonsb_charge
                            contact_count[aicg_itype_sb_qx] += 1
                        else
                            contact_count[aicg_itype_sb_da] += 1
                        end
                    elseif is_nonsb_charge
                        contact_count[aicg_itype_sb_qx] += 1
                    elseif atom_name_1[1] == 'C' || atom_name_2[1] == 'C'
                        contact_count[aicg_itype_sb_cx] += 1
                    else
                        contact_count[aicg_itype_sb_xx] += 1
                    end
                end
            end
        end
    end

    # control the number of long-range contacts
    if aicg_go_atomic_cutoff > aicg_atomic_cutoff
        contact_count[aicg_itype_lr_ct] -= num_short_range_contact
    else
        contact_count[aicg_itype_lr_ct] = 0
    end

    # control the number of salty bridge
    if contact_count[aicg_itype_ss_sb] >= 2
        contact_count[aicg_itype_ss_qx] += contact_count[aicg_itype_ss_sb] - 1
        contact_count[aicg_itype_ss_sb] = 1
    end
    
    return contact_count
end


###############################################################################
#                               Main Function!!!                              #
###############################################################################
function pdb_2_top(pdb_name)

    # ===============
    # Data Structures
    # ===============

    # aa_atom_index  = []
    # aa_atom_name   = []
    # aa_resid_index = []
    # aa_resid_name  = []
    # aa_chain_id    = []
    # aa_coor        = []
    # aa_occupancy   = []
    # aa_tempfactor  = []
    # aa_element     = []
    # aa_charge      = []

    cg_coor            = []
    cg_particle_name   = []
    cg_particle_index  = []
    cg_particle_charge = []
    cg_residue_name    = []
    cg_residue_index   = []

    # ================
    # Step 1: open PDB
    # ================
    println("> Step 1: open PDB file.")

end

function main()
    
    args = parse_commandline()

    pdb_2_top(args["pdb"])

end

main()
