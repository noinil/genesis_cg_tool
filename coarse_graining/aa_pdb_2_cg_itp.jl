#!/usr/bin/env julia

###############################################################################
#                                    README
# This program read PDB structures and prepare toppology and coordinate files
# for CG MD simulations in Genesis.
#
# PDB format:
# 1. Atoms startswith "ATOM  "
# 2. Chains should end with "TER" and have different IDs
###############################################################################

using Printf
using ArgParse
using Formatting
using LinearAlgebra
using ProgressMeter

###########################################################################
#                          Force Field Parameters                         #
###########################################################################
#      ____   _    ____      _    __  __ _____ _____ _____ ____  ____  
#     |  _ \ / \  |  _ \    / \  |  \/  | ____|_   _| ____|  _ \/ ___| 
#     | |_) / _ \ | |_) |  / _ \ | |\/| |  _|   | | |  _| | |_) \___ \ 
#     |  __/ ___ \|  _ <  / ___ \| |  | | |___  | | | |___|  _ < ___) |
#     |_| /_/   \_\_| \_\/_/   \_\_|  |_|_____| |_| |_____|_| \_\____/ 
#     
###########################################################################
 

# ==================
# Physical Constants
# ==================
const CAL2JOU = 4.184

# =====================================
# General Parameters: Mass, Charge, ...
# =====================================

ATOM_MASS_DICT = Dict(
    'C' => 12.011,
    'N' => 14.001,
    'O' => 15.999,
    'P' => 30.974,
    'S' => 32.065,
    'H' =>  1.008
)

RES_MASS_DICT = Dict(
    "ALA" =>  71.09,
    "ARG" => 156.19,
    "ASN" => 114.11,
    "ASP" => 115.09,
    "CYS" => 103.15,
    "CYM" => 103.15,
    "CYT" => 103.15,
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
    "VAL" =>  99.14,
    "DA"  => 134.10,
    "DC"  => 110.10,
    "DG"  => 150.10,
    "DT"  => 125.10,
    "DP"  =>  94.97,
    "DS"  =>  83.11
)

RES_CHARGE_DICT = Dict(
    "ALA" =>  0.0,
    "ARG" =>  1.0,
    "ASN" =>  0.0,
    "ASP" => -1.0,
    "CYS" =>  0.0,
    "CYM" =>  0.0,
    "CYT" =>  0.0,
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
    "VAL" =>  0.0,
    "DA"  =>  0.0,
    "DC"  =>  0.0,
    "DG"  =>  0.0,
    "DT"  =>  0.0,
    "DP"  => -0.6,
    "DS"  =>  0.0
)

RES_NAME_LIST_PROTEIN = (
    "ALA", "ARG", "ASN", "ASP",
    "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    "CYM", "CYT")

RES_NAME_LIST_DNA = ("DA", "DC", "DG", "DT")

RES_NAME_LIST_RNA = ("RA", "RC", "RG", "RU")

# DNA CG residue atom names
ATOM_NAME_LIST_DP = ("P", "OP1", "OP2", "O5'", "O1P", "O2P")
ATOM_NAME_LIST_DS = ("C5'", "C4'", "C3'", "C2'", "C1'", "O4'", "O2'")

# RNA CG residue atom names
ATOM_NAME_LIST_RP = ("P", "OP1", "OP2", "O1P", "O2P")
ATOM_NAME_LIST_RS = ("C5'", "C4'", "C3'", "C2'", "C1'", "O5'", "O4'", "O3'", "O2'")

# ==============
# Molecule Types
# ==============

const MOL_DNA     = 1
const MOL_RNA     = 2
const MOL_PROTEIN = 3
const MOL_OTHER   = 4
MOL_TYPE_LIST = ["DNA", "RNA", "protein", "other", "unknown"]

# ===============================
# Protein AICG2+ Model Parameters
# ===============================

# AICG2+ bond force constant
const AICG_BOND_K               = 110.40 * CAL2JOU * 100 * 2
# AICG2+ sigma for Gaussian angle
const AICG_13_SIGMA             = 0.15 * 0.1  # nm
# AICG2+ sigma for Gaussian dihedral
const AICG_14_SIGMA             = 0.15        # Rad ??
# AICG2+ atomistic contact cutoff
const AICG_GO_ATOMIC_CUTOFF     = 6.5
# AICG2+ pairwise interaction cutoff
const AICG_ATOMIC_CUTOFF        = 5.0
# AICG2+ hydrogen bond cutoff
const AICG_HYDROGEN_BOND_CUTOFF = 3.2
# AICG2+ salt bridge cutoff
const AICG_SALT_BRIDGE_CUTOFF   = 3.5
# AICG2+ energy cutoffs
const AICG_ENE_UPPER_LIM        = -0.5
const AICG_ENE_LOWER_LIM        = -5.0
# average and general AICG2+ energy values
const AICG_13_AVE               = 1.72
const AICG_14_AVE               = 1.23
const AICG_CONTACT_AVE          = 0.55
const AICG_13_GEN               = 1.11
const AICG_14_GEN               = 0.87
const AICG_CONTACT_GEN          = 0.32

# AICG2+ pairwise interaction pairs
const AICG_ITYPE_BB_HB = 1  # B-B hydrogen bonds
const AICG_ITYPE_BB_DA = 2  # B-B donor-accetor contacts
const AICG_ITYPE_BB_CX = 3  # B-B carbon-X contacts
const AICG_ITYPE_BB_XX = 4  # B-B other
const AICG_ITYPE_SS_HB = 5  # S-S hydrogen bonds
const AICG_ITYPE_SS_SB = 6  # S-S salty bridge
const AICG_ITYPE_SS_DA = 7  # S-S donor-accetor contacts
const AICG_ITYPE_SS_CX = 8  # S-S carbon-X contacts
const AICG_ITYPE_SS_QX = 9  # S-S charge-X contacts
const AICG_ITYPE_SS_XX = 10 # S-S other
const AICG_ITYPE_SB_HB = 11 # S-B hydrogen bonds
const AICG_ITYPE_SB_DA = 12 # S-B donor-accetor contacts
const AICG_ITYPE_SB_CX = 13 # S-B carbon-X contacts
const AICG_ITYPE_SB_QX = 14 # S-B charge-X contacts
const AICG_ITYPE_SB_XX = 15 # S-B other
const AICG_ITYPE_LR_CT = 16 # long range contacts
const AICG_ITYPE_OFFST = 17  # offset

AICG_PAIRWISE_ENERGY = zeros(Float64, 17)
AICG_PAIRWISE_ENERGY[AICG_ITYPE_BB_HB] = - 1.4247   # B-B hydrogen bonds
AICG_PAIRWISE_ENERGY[AICG_ITYPE_BB_DA] = - 0.4921   # B-B donor-accetor contacts
AICG_PAIRWISE_ENERGY[AICG_ITYPE_BB_CX] = - 0.2404   # B-B carbon-X contacts
AICG_PAIRWISE_ENERGY[AICG_ITYPE_BB_XX] = - 0.1035   # B-B other
AICG_PAIRWISE_ENERGY[AICG_ITYPE_SS_HB] = - 5.7267   # S-S hydrogen bonds
AICG_PAIRWISE_ENERGY[AICG_ITYPE_SS_SB] = -12.4878   # S-S salty bridge
AICG_PAIRWISE_ENERGY[AICG_ITYPE_SS_DA] = - 0.0308   # S-S donor-accetor contacts
AICG_PAIRWISE_ENERGY[AICG_ITYPE_SS_CX] = - 0.1113   # S-S carbon-X contacts
AICG_PAIRWISE_ENERGY[AICG_ITYPE_SS_QX] = - 0.2168   # S-S charge-X contacts
AICG_PAIRWISE_ENERGY[AICG_ITYPE_SS_XX] =   0.2306   # S-S other
AICG_PAIRWISE_ENERGY[AICG_ITYPE_SB_HB] = - 3.4819   # S-B hydrogen bonds
AICG_PAIRWISE_ENERGY[AICG_ITYPE_SB_DA] = - 0.1809   # S-B donor-accetor contacts
AICG_PAIRWISE_ENERGY[AICG_ITYPE_SB_CX] = - 0.1209   # S-B carbon-X contacts
AICG_PAIRWISE_ENERGY[AICG_ITYPE_SB_QX] = - 0.2984   # S-B charge-X contacts
AICG_PAIRWISE_ENERGY[AICG_ITYPE_SB_XX] = - 0.0487   # S-B other
AICG_PAIRWISE_ENERGY[AICG_ITYPE_LR_CT] = - 0.0395   # long range contacts
AICG_PAIRWISE_ENERGY[AICG_ITYPE_OFFST] = - 0.1051   # offset

# ============================
# DNA 3SPN.2C Model Parameters
# ============================

# 3SPN.2C bond force constant
DNA3SPN_BOND_K_2    = 60.0
DNA3SPN_BOND_K_4    = 600000.0
# 3SPN.2C force constant for Gaussian dihedral
DNA3SPN_DIH_G_K     = 7.0
# 3SPN.2C sigma for Gaussian dihedral
DNA3SPN_DIH_G_SIGMA = 0.3
# 3SPN.2C force constant for Gaussian dihedral
DNA3SPN_DIH_P_K     = 2.0

# ====================
# GRO TOP File Options
# ====================

# "NREXCL" in "[moleculetype]"
const MOL_NR_EXCL             = 3
# "CGNR" in "[atoms]"
const AICG_ATOM_FUNC_NR       = 1
const DNA3SPN_ATOM_FUNC_NR    = 1
# "f" in "[bonds]"
const AICG_BOND_FUNC_TYPE     = 1
const DNA3SPN_BOND_FUNC2_TYPE = 1
const DNA3SPN_BOND_FUNC4_TYPE = 21
# "f" in AICG-type "[angles]"
const AICG_ANG_G_FUNC_TYPE    = 21
# "f" in Flexible-type "[angles]"
const AICG_ANG_F_FUNC_TYPE    = 22
# "f" in DNA "[angles]"
const DNA3SPN_ANG_FUNC_TYPE   = 1
# "f" in AICG-type "[dihedral]"
const AICG_DIH_G_FUNC_TYPE    = 21
# "f" in Flexible-type "[dihedral]"
const AICG_DIH_F_FUNC_TYPE    = 22
# "f" in DNA Gaussian "[dihedral]"
const DNA3SPN_DIH_G_FUNC_TYPE = 21
# "f" in DNA Periodic "[dihedral]"
const DNA3SPN_DIH_P_FUNC_TYPE = 1
const DNA3SPN_DIH_P_FUNC_PERI = 1
# "f" in Go-contacts "[pairs]"
const AICG_CONTACT_FUNC_TYPE  = 2


###############################################################################
#                                  Functions                                  #
###############################################################################
#          ____    _    ____ ___ ____   _____ _   _ _   _  ____ 
#         | __ )  / \  / ___|_ _/ ___| |  ___| | | | \ | |/ ___|
#         |  _ \ / _ \ \___ \| | |     | |_  | | | |  \| | |    
#         | |_) / ___ \ ___) | | |___  |  _| | |_| | |\  | |___ 
#         |____/_/   \_\____/___\____| |_|    \___/|_| \_|\____|
#         
###############################################################################

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


# --------------
# Center of mass
# --------------
function compute_center_of_mass(atom_indices, atom_names, atom_coors)
    total_mass       = 0
    tmp_coor         = zeros(Float64, 3)
    for i in atom_indices
        a_mass       = ATOM_MASS_DICT[atom_names[i][1]]
        a_coor       = atom_coors[:, i]
        total_mass  += a_mass
        tmp_coor    += a_coor
    end
    com = tmp_coor / total_mass
    return com
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
    end
    return false
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
    end
    return false
end
function is_protein_anion(atom_name, res_name)
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
function is_protein_hb_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
    if  is_protein_hb_acceptor(atom_name_1) &&
        is_protein_hb_donor(atom_name_2, res_name_2)
        return true
    elseif is_protein_hb_acceptor(atom_name_2) &&
        is_protein_hb_donor(atom_name_1, res_name_1)
        return true
    end
    return false
end
function is_protein_sb_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
    if  is_protein_cation(atom_name_1, res_name_1) &&
        is_protein_anion(atom_name_2,  res_name_2)
        return true
    elseif is_protein_cation(atom_name_2, res_name_2) &&
        is_protein_anion(atom_name_1,  res_name_1)
        return true
    end
    return false
end
function is_protein_nonsb_charge_pair(atom_name_1, res_name_1, atom_name_2, res_name_2)
    if  is_protein_cation(atom_name_1, res_name_1) ||
        is_protein_anion(atom_name_1,  res_name_1) ||
        is_protein_cation(atom_name_2, res_name_2) ||
        is_protein_anion(atom_name_2,  res_name_2)
        return true
    end
    return false
end
function is_protein_go_contact(resid1, resid2, atom_names, atom_coors)
    for i in resid1.atoms
        atom_name_1 = atom_names[i]
        if atom_name_1[1] == 'H'
            continue
        end
        coor_1 = atom_coors[:, i]
        for j in resid2.atoms
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
function count_aicg_atomic_contact(resid1, resid2, res_name_1, res_name_2, atom_names, atom_coors)
    contact_count                   = zeros(Int, 17)
    contact_count[AICG_ITYPE_OFFST] = 1
    num_short_range_contact         = 0
    for i in resid1.atoms
        atom_name_1 = atom_names[i]
        if atom_name_1[1] == 'H'
            continue
        end
        coor_1 = atom_coors[:, i]
        for j in resid2.atoms
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
        contact_count[AICG_ITYPE_LR_CT] = 0
    end

    # control the number of salty bridge
    if contact_count[AICG_ITYPE_SS_SB] >= 2
        contact_count[AICG_ITYPE_SS_QX] += contact_count[AICG_ITYPE_SS_SB] - 1
        contact_count[AICG_ITYPE_SS_SB] = 1
    end

    return contact_count
end

# -----------------
# 3SPN.2C DNA model
# -----------------



###############################################################################
#                                Core Function                                #
###############################################################################
struct AAResidue
    name::String
    atoms::Array{Int64, 1}
end

struct AAChain
    id::Char
    residues::Array{Int64, 1}
end

struct CGResidue
    res_idx::Int
    res_name::String
    atm_name::String
    atoms::Array{Int64, 1}
end

struct CGChain
    first::Int
    last::Int
    moltype::Int
end

###############################################################################
#                            ____ ___  ____  _____ 
#                           / ___/ _ \|  _ \| ____|
#                          | |  | | | | |_) |  _|  
#                          | |__| |_| |  _ <| |___ 
#                           \____\___/|_| \_\_____|
# 
###############################################################################
function pdb_2_top(pdb_name, protein_charge_filename, scale_scheme)

    # ===============
    # Data Structures
    # ===============

    aa_num_atom = 0
    aa_num_residue = 0
    aa_num_chain = 0

    # ================
    # Step 1: open PDB
    # ================
    println("============================================================")
    println("> Step 1: open PDB file.")

    aa_pdb_lines = []
    for line in eachline(pdb_name)
        if startswith(line, "ATOM")
            push!(aa_pdb_lines, rpad(line, 80))
            aa_num_atom += 1
        elseif startswith(line, "TER")
            push!(aa_pdb_lines, rpad(line, 80))
        end
    end

    aa_atom_name   = fill("    ",       aa_num_atom)
    aa_coor        = zeros(Float64, (3, aa_num_atom))

    aa_residues = []
    aa_chains   = []

    i_atom     = 0
    i_resid    = 0
    curr_resid = NaN
    curr_chain = NaN
    curr_rname = "    "
    residue_name = "    "
    chain_id = '?'
    tmp_res_atoms = []
    tmp_chain_res = []
    for line in aa_pdb_lines
        if startswith(line, "TER")
            if length(tmp_res_atoms) > 0
                push!(aa_residues, AAResidue(residue_name, tmp_res_atoms))
                tmp_res_atoms = []
            end
            if length(tmp_chain_res) > 0
                push!(aa_chains, AAChain(chain_id, tmp_chain_res))
                tmp_chain_res = []
            end
            continue
        end

        i_atom += 1
        atom_serial       = parse(Int, line[7:11])
        atom_name         = strip(line[13:16])
        residue_name      = strip(line[18:21])
        chain_id          = line[22]
        residue_serial    = parse(Int, line[23:26])
        coor_x            = parse(Float64, line[31:38])
        coor_y            = parse(Float64, line[39:46])
        coor_z            = parse(Float64, line[47:54])
        # occupancy         = parse(Float64, line[55:60])
        # tempfactor        = parse(Float64, line[61:66])
        # segment_id        = strip(line[67:76])
        # element_name      = strip(line[77:78])
        # charge            = parse(Float64, line[79:80])

        aa_atom_name[i_atom]   = atom_name
        aa_coor[1, i_atom]     = coor_x
        aa_coor[2, i_atom]     = coor_y
        aa_coor[3, i_atom]     = coor_z

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
    aa_num_residue = length(aa_residues)
    aa_num_chain = length(aa_chains)
    println("          > Number of atoms:    $(aa_num_atom)")
    println("          > Number of residues: $(aa_num_residue)")
    println("          > Number of chains:   $(aa_num_chain)")

    # ===============================
    # Step 2: find out molecule types
    # ===============================
    println("============================================================")
    println("> Step 2: set molecular types for every chain.")

    cg_num_particles = 0

    cg_chain_mol_types = zeros(Int, aa_num_chain)
    cg_chain_length    = zeros(Int, aa_num_chain)

    for i_chain = 1:aa_num_chain
        chain = aa_chains[i_chain]
        mol_type = -1
        for i_res in chain.residues
            res_name = aa_residues[i_res].name
            if in(res_name, RES_NAME_LIST_PROTEIN)
                tmp_mol_type = MOL_PROTEIN
            elseif in(res_name, RES_NAME_LIST_DNA)
                tmp_mol_type = MOL_DNA
            elseif in(res_name, RES_NAME_LIST_RNA)
                tmp_mol_type = MOL_RNA
            else
                tmp_mol_type = MOL_OTHER
            end
            if mol_type == -1
                mol_type = tmp_mol_type
            elseif tmp_mol_type != mol_type
                error("BUG: Inconsistent residue types in chain ", i_chain,
                      " ID - ", chain.id,
                      " residue - ", i_res,
                      " : ", res_name)
            end
        end
        cg_chain_mol_types[i_chain] = mol_type
        n_res = length(chain.residues)
        if mol_type == MOL_DNA
            n_particles = 3 * n_res - 1
        elseif mol_type == MOL_RNA
            n_particles = 3 * n_res - 1
        elseif mol_type == MOL_PROTEIN
            n_particles = n_res
        else
            n_particles = 0
        end
        cg_chain_length[i_chain] = n_particles
        cg_num_particles += n_particles
        @printf("          > Chain %3d | %7s \n",
                i_chain, MOL_TYPE_LIST[ mol_type ])
    end

    # ===========================
    # Step 3: Assign CG particles
    # ===========================
    println("============================================================")
    println("> Step 3: assign coarse-grained particles.")

    cg_residues = []
    cg_chains   = []

    i_offset_cg_particle = 0
    i_offset_cg_residue  = 0
    for i_chain in 1:aa_num_chain
        chain = aa_chains[i_chain]
        mol_type = cg_chain_mol_types[i_chain]

        i_bead = i_offset_cg_particle
        i_resi = i_offset_cg_residue
        if mol_type == MOL_PROTEIN
            for i_res in chain.residues
                cg_idx = []
                res_name = aa_residues[i_res].name
                for i_atom in aa_residues[i_res].atoms
                    atom_name = aa_atom_name[i_atom]
                    if atom_name[1] == 'H'
                        continue
                    else
                        push!(cg_idx, i_atom)
                    end
                end
                i_bead += 1
                i_resi += 1
                push!(cg_residues, CGResidue(i_resi, res_name, "CA", cg_idx))
            end
        elseif mol_type == MOL_DNA
            tmp_atom_index_O3p = 0
            for (i_local_index, i_res) in enumerate( chain.residues )
                res_name = aa_residues[i_res].name
                cg_DP_idx = [tmp_atom_index_O3p]
                cg_DS_idx = []
                cg_DB_idx = []
                for i_atom in aa_residues[i_res].atoms
                    atom_name = aa_atom_name[i_atom]
                    if atom_name[1] == 'H'
                        continue
                    elseif in(atom_name, ATOM_NAME_LIST_DP)
                        push!(cg_DP_idx, i_atom)
                    elseif in(atom_name, ATOM_NAME_LIST_DS)
                        push!(cg_DS_idx, i_atom)
                    elseif atom_name == "O3'"
                        tmp_atom_index_O3p = i_atom
                    else
                        push!(cg_DB_idx, i_atom)
                    end
                end
                i_resi += 1
                if i_local_index > 1
                    i_bead += 1
                    push!(cg_residues, CGResidue(i_resi, res_name, "DP", cg_DP_idx))
                end
                i_bead += 1
                push!(cg_residues, CGResidue(i_resi, res_name, "DS", cg_DS_idx))
                i_bead += 1
                push!(cg_residues, CGResidue(i_resi, res_name, "DB", cg_DB_idx))
            end
        elseif mol_type == MOL_RNA
            for (i_local_index, i_res) in enumerate( chain.residues )
                res_name = aa_residues[i_res].name
                cg_RP_idx = []
                cg_RS_idx = []
                cg_RB_idx = []
                for i_atom in aa_residues[i_res].atoms
                    atom_name = aa_atom_name[i_atom]
                    if atom_name[1] == 'H'
                        continue
                    elseif in(atom_name, ATOM_NAME_LIST_RP)
                        push!(cg_RP_idx, i_atom)
                    elseif in(atom_name, ATOM_NAME_LIST_RS)
                        push!(cg_RS_idx, i_atom)
                    else
                        push!(cg_RB_idx, i_atom)
                    end
                end
                i_resi += 1
                if i_local_index > 1
                    i_bead += 1
                    push!(cg_residues, CGResidue(i_resi, res_name, "RP", cg_RP_idx))
                end
                i_bead += 1
                push!(cg_residues, CGResidue(i_resi, res_name, "RS", cg_RS_idx))
                i_bead += 1
                push!(cg_residues, CGResidue(i_resi, res_name, "RB", cg_RB_idx))
            end
        end
        push!(cg_chains, CGChain(i_offset_cg_particle + 1, i_bead, mol_type))
        i_offset_cg_particle += cg_chain_length[i_chain]
        i_offset_cg_residue  += length(chain.residues)
    end
    for i_chain in 1:aa_num_chain
        @printf("          > Chain %3d | # particles: %5d | %5d -- %5d \n",
                i_chain, cg_chain_length[i_chain],
                cg_chains[i_chain].first, cg_chains[i_chain].last)
    end
    println("------------------------------------------------------------")
    println("          In total: $(cg_num_particles) CG particles.")




    # =========================================================================
    #        ____ ____   _____ ___  ____   ___  _     ___   ______   __
    #       / ___/ ___| |_   _/ _ \|  _ \ / _ \| |   / _ \ / ___\ \ / /
    #      | |  | |  _    | || | | | |_) | | | | |  | | | | |  _ \ V / 
    #      | |__| |_| |   | || |_| |  __/| |_| | |__| |_| | |_| | | |  
    #       \____\____|   |_| \___/|_|    \___/|_____\___/ \____| |_|  
    # 
    # =========================================================================

    cg_resid_name  = fill("    ", cg_num_particles)
    cg_resid_index = zeros(Int, cg_num_particles)
    cg_bead_name   = fill("    ", cg_num_particles)
    cg_bead_charge = zeros(Float64, cg_num_particles)
    cg_bead_mass   = zeros(Float64, cg_num_particles)
    cg_bead_coor   = zeros(Float64, (3, cg_num_particles))

    # protein
    top_cg_pro_bonds        = []
    top_cg_pro_angles       = []
    top_cg_pro_dihedrals    = []
    top_cg_pro_aicg13       = []
    top_cg_pro_aicg14       = []
    top_cg_pro_aicg_contact = []

    param_cg_pro_e_13       = []
    param_cg_pro_e_14       = []
    param_cg_pro_e_contact  = []

    # DNA
    top_cg_DNA_bonds        = []
    top_cg_DNA_angles       = []
    top_cg_DNA_dih_Gaussian = []
    top_cg_DNA_dih_periodic = []

    # =================================
    # Step 4: AICG2+ model for proteins
    # =================================
    #                  _       _       
    #  _ __  _ __ ___ | |_ ___(_)_ __  
    # | '_ \| '__/ _ \| __/ _ \ | '_ \ 
    # | |_) | | | (_) | ||  __/ | | | |
    # | .__/|_|  \___/ \__\___|_|_| |_|
    # |_|  
    # 
    # =================================
    println("============================================================")
    println("> Step 4: processing proteins.")

    # --------------------------------
    # Step 4.1: find out C-alpha atoms
    # --------------------------------
    println("------------------------------------------------------------")
    println(">      4.1: determine CA mass, charge, and coordinates.")

    for i_chain in 1:aa_num_chain
        chain = cg_chains[i_chain]

        if chain.moltype != MOL_PROTEIN
            continue
        end

        for i_res in chain.first : chain.last
            res_name = cg_residues[i_res].res_name
            for i_atom in cg_residues[i_res].atoms
                if aa_atom_name[i_atom] == "CA"
                    cg_resid_name[i_res]   = res_name
                    cg_resid_index[i_res]  = cg_residues[i_res].res_idx
                    cg_bead_name[i_res]    = "CA"
                    cg_bead_charge[i_res]  = RES_CHARGE_DICT[res_name]
                    cg_bead_mass[i_res]    = RES_MASS_DICT[res_name]
                    cg_bead_coor[:, i_res] = aa_coor[:, i_atom]
                    break
                end
            end
        end
    end

    if length(protein_charge_filename) > 0
        try
            for line in eachline(protein_charge_filename)
                charge_data = split(line)
                if length(charge_data) < 1
                    continue
                end
                i = parse(Int, charge_data[1])
                c = parse(Float64, charge_data[2])
                cg_bead_charge[i] = c
            end
        catch e
            println(e)
            error("ERROR in user-defined charge distribution.\n")
        end
    end
    println(">           ... DONE!")

    # -------------------------
    # Step 4.2: AICG2+ topology
    # -------------------------
    println("------------------------------------------------------------")
    println(">      4.2: AICG2+ topology.")
    println(" - - - - - - - - - - - - - - - - - - - - - - - -")
    println(">      4.2.1: AICG2+ local interactions.")
    for i_chain in 1:aa_num_chain
        chain = cg_chains[i_chain]

        if chain.moltype != MOL_PROTEIN
            continue
        end

        for i_res in chain.first : chain.last - 1
            coor1 = cg_bead_coor[:, i_res]
            coor2 = cg_bead_coor[:, i_res + 1]
            dist12 = compute_distance(coor1, coor2)
            push!(top_cg_pro_bonds, (i_res, dist12))
        end
    end
    println(">           ... Bond: DONE!")

    e_ground_local = 0.0
    e_ground_13    = 0.0
    num_angle      = 0
    for i_chain in 1:aa_num_chain
        chain = cg_chains[i_chain]

        if chain.moltype != MOL_PROTEIN
            continue
        end

        for i_res in chain.first : chain.last - 2
            coor1 = cg_bead_coor[:, i_res]
            coor3 = cg_bead_coor[:, i_res + 2]
            dist13 = compute_distance(coor1, coor3)
            push!(top_cg_pro_angles, i_res)
            push!(top_cg_pro_aicg13, (i_res, dist13))

            # count AICG2+ atomic contact
            contact_counts = count_aicg_atomic_contact(cg_residues[ i_res ],
                                                       cg_residues[ i_res + 2 ],
                                                       cg_resid_name[i_res],
                                                       cg_resid_name[i_res + 2],
                                                       aa_atom_name,
                                                       aa_coor)

            # calculate AICG2+ pairwise energy
            e_local = dot(AICG_PAIRWISE_ENERGY, contact_counts)
            if e_local > AICG_ENE_UPPER_LIM
                e_local = AICG_ENE_UPPER_LIM
            end
            if e_local < AICG_ENE_LOWER_LIM
                e_local = AICG_ENE_LOWER_LIM
            end
            e_ground_local += e_local
            e_ground_13    += e_local
            num_angle      += 1
            push!(param_cg_pro_e_13, e_local)
        end
    end
    println(">           ... Angle: DONE!")

    e_ground_14 = 0.0
    num_dih = 0
    for i_chain in 1:aa_num_chain
        chain = cg_chains[i_chain]

        if chain.moltype != MOL_PROTEIN
            continue
        end

        for i_res in chain.first : chain.last - 3
            coor1   = cg_bead_coor[:, i_res]
            coor2   = cg_bead_coor[:, i_res + 1]
            coor3   = cg_bead_coor[:, i_res + 2]
            coor4   = cg_bead_coor[:, i_res + 3]
            dihed = compute_dihedral(coor1, coor2, coor3, coor4)
            push!(top_cg_pro_dihedrals, i_res)
            push!(top_cg_pro_aicg14, (i_res, dihed))

            # count AICG2+ atomic contact
            contact_counts = count_aicg_atomic_contact(cg_residues[ i_res ],
                                                       cg_residues[ i_res + 3 ],
                                                       cg_resid_name[i_res],
                                                       cg_resid_name[i_res + 3],
                                                       aa_atom_name,
                                                       aa_coor)

            # calculate AICG2+ pairwise energy
            e_local = dot(AICG_PAIRWISE_ENERGY, contact_counts)
            if e_local > AICG_ENE_UPPER_LIM
                e_local = AICG_ENE_UPPER_LIM
            end
            if e_local < AICG_ENE_LOWER_LIM
                e_local = AICG_ENE_LOWER_LIM
            end
            e_ground_local += e_local
            e_ground_14    += e_local
            num_dih      += 1
            push!(param_cg_pro_e_14, e_local)
        end
    end
    println(">           ... Dihedral: DONE!")

    # ------------------------
    # Normalize local energies
    # ------------------------
    e_ground_local /= (num_angle + num_dih)
    e_ground_13    /= num_angle
    e_ground_14    /= num_dih

    if scale_scheme == 0
        for i in 1:length(param_cg_pro_e_13)
            param_cg_pro_e_13[i] *= AICG_13_AVE / e_ground_13
        end
        for i in 1:length(param_cg_pro_e_14)
            param_cg_pro_e_14[i] *= AICG_14_AVE / e_ground_14
        end
    elseif scale_scheme == 1
        for i in 1:length(param_cg_pro_e_13)
            param_cg_pro_e_13[i] *= -AICG_13_GEN
        end
        for i in 1:length(param_cg_pro_e_14)
            param_cg_pro_e_14[i] *= -AICG_14_GEN
        end
    end



    # -----------------------
    # Go type native contacts
    # -----------------------
    println(" - - - - - - - - - - - - - - - - - - - - - - - -")
    println(">      4.2.2: AICG2+ Go-type native contacts.")
    e_ground_contact = 0.0
    num_contact = 0

    # intra-molecular contacts
    @showprogress 1 "        Calculating intra-molecular contacts..."   for i_chain in 1:aa_num_chain
        chain = cg_chains[i_chain]

        if chain.moltype != MOL_PROTEIN
            continue
        end

        for i_res in chain.first : chain.last - 4
            coor_cai = cg_bead_coor[:, i_res]
            for j_res in i_res + 4 : chain.last
                coor_caj = cg_bead_coor[:, j_res]
                if is_protein_go_contact(cg_residues[i_res], cg_residues[j_res], aa_atom_name, aa_coor)
                    native_dist = compute_distance(coor_cai, coor_caj)
                    num_contact += 1
                    push!(top_cg_pro_aicg_contact, (i_res, j_res, native_dist))

                    # count AICG2+ atomic contact
                    contact_counts = count_aicg_atomic_contact(cg_residues[ i_res ],
                                                               cg_residues[ j_res ],
                                                               cg_resid_name[i_res],
                                                               cg_resid_name[j_res],
                                                               aa_atom_name,
                                                               aa_coor)

                    # calculate AICG2+ pairwise energy
                    e_local = dot(AICG_PAIRWISE_ENERGY, contact_counts)
                    if e_local > AICG_ENE_UPPER_LIM
                        e_local = AICG_ENE_UPPER_LIM
                    end
                    if e_local < AICG_ENE_LOWER_LIM
                        e_local = AICG_ENE_LOWER_LIM
                    end
                    e_ground_contact += e_local
                    num_contact      += 1
                    push!(param_cg_pro_e_contact, e_local)
                end
            end
        end
    end
    println(">           ... intra-molecular contacts: DONE!")

    # inter-molecular ( protein-protein ) contacts
    @showprogress 1 "        Calculating inter-molecular contacts..." for i_chain in 1 : aa_num_chain - 1
        chain1 = cg_chains[i_chain]

        if chain1.moltype != MOL_PROTEIN
            continue
        end

        for j_chain in i_chain + 1 : aa_num_chain
            chain2 = cg_chains[j_chain]

            if chain2.moltype != MOL_PROTEIN
                continue
            end

            for i_res in chain1.first : chain1.last
                coor_cai = cg_bead_coor[:, i_res]
                for j_res in chain2.first : chain2.last
                    coor_caj = cg_bead_coor[:, j_res]
                    if is_protein_go_contact(cg_residues[i_res], cg_residues[j_res], aa_atom_name, aa_coor)
                        native_dist = compute_distance(coor_cai, coor_caj)
                        num_contact += 1
                        push!(top_cg_pro_aicg_contact, (i_res, j_res, native_dist))

                        # count AICG2+ atomic contact
                        contact_counts = count_aicg_atomic_contact(cg_residues[ i_res ],
                                                                   cg_residues[ j_res ],
                                                                   cg_resid_name[i_res],
                                                                   cg_resid_name[j_res],
                                                                   aa_atom_name,
                                                                   aa_coor)

                        # calculate AICG2+ pairwise energy
                        e_local = dot(AICG_PAIRWISE_ENERGY, contact_counts)
                        if e_local > AICG_ENE_UPPER_LIM
                            e_local = AICG_ENE_UPPER_LIM
                        end
                        if e_local < AICG_ENE_LOWER_LIM
                            e_local = AICG_ENE_LOWER_LIM
                        end
                        e_ground_contact += e_local
                        num_contact      += 1
                        push!(param_cg_pro_e_contact, e_local)
                    end
                end
            end
        end
    end
    println(">           ... inter-molecular contacts: DONE!")

    # normalize
    e_ground_contact /= num_contact

    if scale_scheme == 0
        for i in 1:length(param_cg_pro_e_contact)
            param_cg_pro_e_contact[i] *= AICG_CONTACT_AVE / e_ground_contact
        end
    elseif scale_scheme == 1
        for i in 1:length(param_cg_pro_e_contact)
            param_cg_pro_e_contact[i] *= -AICG_CONTACT_GEN
        end
    end

    println("------------------------------------------------------------")
    @printf("          > Total number of protein contacts: %12d  \n",
            length( top_cg_pro_aicg_contact ))





    # =============================
    # Step 5: 3SPN.2C model for DNA
    # =============================
    #         _             
    #      __| |_ __   __ _ 
    #     / _` | '_ \ / _` |
    #    | (_| | | | | (_| |
    #     \__,_|_| |_|\__,_|
    #    
    # =============================
    println("============================================================")
    println("> Step 5: processing DNA.")

    # ---------------------------
    # Step 5.1: determine P, S, B
    # ---------------------------
    println("------------------------------------------------------------")
    println(">      5.1: determine P, S, B mass, charge, and coordinates.")

    for i_chain in 1 : aa_num_chain
        chain = cg_chains[i_chain]

        if chain.moltype != MOL_DNA
            continue
        end

        for i_res in chain.first : chain.last
            res_name  = cg_residues[i_res].res_name
            bead_name = cg_residues[i_res].atm_name
            bead_coor = compute_center_of_mass(cg_residues[i_res].atoms, aa_atom_name, aa_coor)
            cg_resid_name[i_res]   = res_name
            cg_resid_index[i_res]  = cg_residues[i_res].res_idx
            cg_bead_name[i_res]    = bead_name
            cg_bead_charge[i_res]  = RES_CHARGE_DICT[res_name]
            cg_bead_mass[i_res]    = RES_MASS_DICT[res_name]
            cg_bead_coor[:, i_res] = bead_coor
        end
    end

    println(">           ... DONE!")







    # =======================
    # Output: itp & gro files
    # =======================
    println("============================================================")
    println("> Step 8: output .itp and .gro files.")

    # ------------------------------------------------------------
    # itp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ------------------------------------------------------------
    itp_mol_head = "[ moleculetype ]\n"
    itp_mol_comm = format(";{1:15s} {2:6s}\n", "name", "nrexcl")
    itp_mol_line = "{1:<16} {2:>6d}\n"

    itp_atm_head = "[ atoms ]\n"
    itp_atm_comm = format(";{:>9}{:>5}{:>10}{:>5}{:>5}{:>5} {:>8} {:>8}\n", "nr", "type", "resnr", "res", "atom", "cg", "charge", "mass")
    itp_atm_line = "{:>10d}{:>5}{:>10d}{:>5}{:>5}{:>5d} {:>8.3f} {:>8.3f}\n"

    itp_bnd_head = "[ bonds ]\n"
    itp_bnd_comm = format(";{:>9}{:>10}{:>5}{:>18}{:>18}\n", "i", "j", "f", "eq", "k2")
    itp_bnd_line = "{:>10d}{:>10d}{:>5d}{:>18.4E}{:>18.4E}\n"

    itp_13_head = "[ angles ] ; AICG2+ 1-3 interaction\n"
    itp_13_comm = format(";{:>9}{:>10}{:>10}{:>5}{:>15}{:>15}{:>15}\n", "i", "j", "k", "f", "eq", "k", "w")
    itp_13_line = "{:>10d}{:>10d}{:>10d}{:>5d}{:>15.4E}{:>15.4E}{:>15.4E}\n"

    itp_ang_head = "[ angles ] ; AICG2+ flexible local interaction\n"
    itp_ang_comm = format(";{:>9}{:>10}{:>10}{:>5}\n", "i", "j", "k", "f")
    itp_ang_line = "{:>10d}{:>10d}{:>10d}{:>5d}\n"

    itp_dih_G_head = "[ dihedrals ] ; AICG2+ Gaussian dihedrals\n"
    itp_dih_G_comm = format(";{:>9}{:>10}{:>10}{:>10}{:>5}{:>15}{:>15}{:>15}\n", "i", "j", "k", "l", "f", "eq", "k", "w")
    itp_dih_G_line = "{:>10d}{:>10d}{:>10d}{:>10d}{:>5d}{:>15.4E}{:>15.4E}{:>15.4E}\n"

    itp_dih_F_head = "[ dihedrals ] ; AICG2+ flexible local interation\n"
    itp_dih_F_comm = format(";{:>9}{:>10}{:>10}{:>10}{:>5}\n", "i", "j", "k", "l", "f")
    itp_dih_F_line = "{:>10d}{:>10d}{:>10d}{:>10d}{:>5d}\n"

    itp_contact_head = "[ pairs ] ; Go-type native contact\n"
    itp_contact_comm = format(";{:>9}{:>10}{:>10}{:>15}{:>15}\n", "i", "j", "f", "eq", "k")
    itp_contact_line = "{:>10d}{:>10d}{:>10d}{:>15.4E}{:>15.4E}\n"

    itp_exc_head = "[ exclusions ] ; Genesis exclusion list\n"
    itp_exc_comm = format(";{:>9}{:>10}\n", "i", "j")
    itp_exc_line = "{:>10d}{:>10d}\n"

    # --------
    # filename
    # --------
    itp_name = pdb_name[1:end-4] * ".itp"
    itp_file = open(itp_name, "w")

    # --------------------
    # Writing CG particles
    # --------------------
    # write molecule type information
    itp_system_name = pdb_name[1:end-4]
    write(itp_file, itp_mol_head)
    write(itp_file, itp_mol_comm)
    printfmt(itp_file, itp_mol_line, itp_system_name, MOL_NR_EXCL)
    write(itp_file,"\n")

    # write atoms information
    write(itp_file, itp_atm_head)
    write(itp_file, itp_atm_comm)
    for i_bead in 1 : cg_num_particles
        printfmt(itp_file,
                 itp_atm_line,
                 i_bead,
                 cg_resid_name[i_bead],
                 cg_resid_index[i_bead],
                 cg_resid_name[i_bead],
                 cg_bead_name[i_bead],
                 AICG_ATOM_FUNC_NR,
                 cg_bead_charge[i_bead],
                 cg_bead_mass[i_bead])
    end
    write(itp_file,"\n")

    # write bond information
    write(itp_file, itp_bnd_head)
    write(itp_file, itp_bnd_comm)
    for i_bond in 1 : length(top_cg_pro_bonds)
        printfmt(itp_file,
                 itp_bnd_line,
                 top_cg_pro_bonds[i_bond][1],
                 top_cg_pro_bonds[i_bond][1] + 1,
                 AICG_BOND_FUNC_TYPE,
                 top_cg_pro_bonds[i_bond][2] * 0.1,
                 AICG_BOND_K)
    end
    write(itp_file, "\n")

    # write 13 interaction information
    write(itp_file, itp_13_head)
    write(itp_file, itp_13_comm)
    for i_13 in 1 : length(top_cg_pro_aicg13)
        printfmt(itp_file,
                 itp_13_line,
                 top_cg_pro_aicg13[i_13][1],
                 top_cg_pro_aicg13[i_13][1] + 1,
                 top_cg_pro_aicg13[i_13][1] + 2,
                 AICG_ANG_G_FUNC_TYPE,
                 top_cg_pro_aicg13[i_13][2] * 0.1,
                 param_cg_pro_e_13[i_13] * CAL2JOU,
                 AICG_13_SIGMA)
    end
    write(itp_file, "\n")

    # write angle interaction information
    write(itp_file, itp_ang_head)
    write(itp_file, itp_ang_comm)
    for i_ang in 1 : length(top_cg_pro_angles)
        printfmt(itp_file,
                 itp_ang_line,
                 top_cg_pro_angles[i_ang],
                 top_cg_pro_angles[i_ang] + 1,
                 top_cg_pro_angles[i_ang] + 2,
                 AICG_ANG_F_FUNC_TYPE)
    end
    write(itp_file, "\n")

    # write Gaussian dihedral information
    write(itp_file, itp_dih_G_head)
    write(itp_file, itp_dih_G_comm)
    for i_dih in 1 : length(top_cg_pro_aicg14)
        printfmt(itp_file,
                 itp_dih_G_line,
                 top_cg_pro_aicg14[i_dih][1],
                 top_cg_pro_aicg14[i_dih][1] + 1,
                 top_cg_pro_aicg14[i_dih][1] + 2,
                 top_cg_pro_aicg14[i_dih][1] + 3,
                 AICG_DIH_G_FUNC_TYPE,
                 top_cg_pro_aicg14[i_dih][2],
                 param_cg_pro_e_14[i_dih] * CAL2JOU,
                 AICG_14_SIGMA)
    end
    write(itp_file, "\n")

    # write local flexible dihedral information
    write(itp_file, itp_dih_F_head)
    write(itp_file, itp_dih_F_comm)
    for i_dih in 1 : length(top_cg_pro_dihedrals)
        printfmt(itp_file,
                 itp_dih_F_line,
                 top_cg_pro_dihedrals[i_dih],
                 top_cg_pro_dihedrals[i_dih] + 1,
                 top_cg_pro_dihedrals[i_dih] + 2,
                 top_cg_pro_dihedrals[i_dih] + 3,
                 AICG_DIH_F_FUNC_TYPE)
    end
    write(itp_file, "\n")

    # write Go-type native contacts
    write(itp_file, itp_contact_head)
    write(itp_file, itp_contact_comm)
    for i_c in 1 : length(top_cg_pro_aicg_contact)
        printfmt(itp_file,
                 itp_contact_line,
                 top_cg_pro_aicg_contact[i_c][1],
                 top_cg_pro_aicg_contact[i_c][2],
                 AICG_CONTACT_FUNC_TYPE,
                 top_cg_pro_aicg_contact[i_c][3] * 0.1,
                 param_cg_pro_e_contact[i_c] * CAL2JOU)
    end
    write(itp_file, "\n")

    # write Genesis local-exclusion list
    write(itp_file, itp_exc_head)
    write(itp_file, itp_exc_comm)
    for i_c in 1 : length(top_cg_pro_aicg_contact)
        printfmt(itp_file,
                 itp_exc_line,
                 top_cg_pro_aicg_contact[i_c][1],
                 top_cg_pro_aicg_contact[i_c][2])
    end
    write(itp_file, "\n")

    close(itp_file)
    println(">           ... .itp: DONE!")

    # ------------------------------------------------------------
    # itp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ------------------------------------------------------------
    # HEAD: time in the unit of ps
    GRO_HEAD_STR  = "{}, t= {:>16.3f} \n"
    # ATOM NUM: free format int
    GRO_ATOM_NUM  = "{:>12d} \n"
    # XYZ: in the unit of nm!!!
    GRO_ATOM_LINE = "{:>5d}{:>5}{:>5}{:>5d} {:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f} \n"
    GRO_BOX_LINE  = "{:>15.4f}{:>15.4f}{:>15.4f} \n\n"

    gro_name = pdb_name[1:end-4] * ".gro"
    gro_file = open(gro_name, "w")

    printfmt(gro_file, GRO_HEAD_STR, "CG model for GENESIS: ", 0)
    printfmt(gro_file, GRO_ATOM_NUM, cg_num_particles)

    for i_bead in 1 : cg_num_particles
        printfmt(gro_file, GRO_ATOM_LINE,
                 cg_resid_index[i_bead],
                 cg_resid_name[i_bead],
                 cg_bead_name[i_bead],
                 i_bead,
                 cg_bead_coor[1 , i_bead] * 0.1,
                 cg_bead_coor[2 , i_bead] * 0.1,
                 cg_bead_coor[3 , i_bead] * 0.1,
                 0.0, 0.0, 0.0)
    end
    printfmt(gro_file, GRO_BOX_LINE, 0.0, 0.0, 0.0)

    close(gro_file)
    println(">           ... .gro: DONE!")
    println("------------------------------------------------------------")

    println("[1;32m DONE! [0m ")
    println(" Please check the .itp and .gro files.")
    println("============================================================")
end

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

        "--respac", "-c"
        help = "RESPAC protein charge distribution data."
        arg_type = String
        default = ""

        "--aicg-scale"
        help = "Scale AICG2+ local interactions: 0) average; 1) general (default)."
        arg_type = Int
        default = 1
    end

    return parse_args(s)
end

# ====
# Main
# ====

function main()

    args = parse_commandline()

    pdb_2_top(args["pdb"], args["respac"], args["aicg-scale"])

end

main()
