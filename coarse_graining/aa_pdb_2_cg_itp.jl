#!/usr/bin/env julia

###############################################################################
#                                    README
# 
# This program read PDB structures and prepare toppology and coordinate files
# for CG MD simulations in Genesis.
#
# PDB format:
# 1. Atoms startswith "ATOM  "
# 2. Chains should end with "TER" and have different IDs
# 
# Unit in the script: kcal/mol
# Unit for output:    kJ/mol
###############################################################################

using Printf
using ArgParse
using Formatting
using LinearAlgebra

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
const JOU2CAL = 1.0 / CAL2JOU

# ============
# Force fields
# ============

# protein
const FF_pro_AICG2p      = 1
const FF_pro_Clementi_Go = 2
const FF_pro_KB_Go       = 3
# DNA
const FF_DNA_3SPN2C      = 1
# RNA
const FF_RNA_Go          = 1
# unknown
const FF_UNKNOWN         = 0

FF_PRO_DICT = Dict(
    "AICG2+"   => FF_pro_AICG2p,
    "Clementi" => FF_pro_Clementi_Go,
    "KB-Go"    => FF_pro_KB_Go
)

FF_DNA_DICT = Dict(
    "3SPN.2C"  => FF_DNA_3SPN2C
)

FF_RNA_DICT = Dict(
    "Go"       => FF_RNA_Go
)

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
    "UNK" => 100.00,
    "DA"  => 134.10,
    "DC"  => 110.10,
    "DG"  => 150.10,
    "DT"  => 125.10,
    "DP"  =>  94.97,
    "DS"  =>  83.11,
    "RA"  => 134.10,
    "RC"  => 110.10,
    "RG"  => 150.10,
    "RU"  => 111.10,
    "RP"  =>  62.97,
    "RS"  => 131.11
)

RES_CHARGE_DICT = Dict(
    "ALA" =>  0.0,
    "ARG" =>  1.0,
    "ASN" =>  0.0,
    "ASP" => -1.0,
    "CYS" =>  0.0,
    "CYM" =>  0.0,
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
    "UNK" =>  0.0,
    "DA"  =>  0.0,
    "DC"  =>  0.0,
    "DG"  =>  0.0,
    "DT"  =>  0.0,
    "DP"  => -0.6,
    "DS"  =>  0.0,
    "RA"  =>  0.0,
    "RC"  =>  0.0,
    "RG"  =>  0.0,
    "RU"  =>  0.0,
    "RP"  => -1.0,
    "RS"  =>  0.0
)

RES_SHORTNAME_DICT = Dict(
    "ALA" => "A",
    "ARG" => "R",
    "ASN" => "N",
    "ASP" => "D",
    "CYS" => "C",
    "CYM" => "C",
    "GLN" => "Q",
    "GLU" => "E",
    "GLY" => "G",
    "HIS" => "H",
    "ILE" => "I",
    "LEU" => "L",
    "LYS" => "K",
    "MET" => "M",
    "PHE" => "F",
    "PRO" => "P",
    "SER" => "S",
    "THR" => "T",
    "TRP" => "W",
    "TYR" => "Y",
    "VAL" => "V",
    "UNK" => "X",
    "DA"  => "A",
    "DC"  => "C",
    "DG"  => "G",
    "DT"  => "T",
    "RA"  => "A",
    "RC"  => "C",
    "RG"  => "G",
    "RU"  => "U"
)

RES_NAME_LIST_PROTEIN = (
    "ALA", "ARG", "ASN", "ASP",
    "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    "CYM", "UNK")

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

# ====================================
# Protein Clementi Go Model Parameters
# ====================================

# Clementi Go energy unit: epsilon
const CCGO_EPSILON           = 1.0
# Clementi Go bond force constant
const CCGO_BOND_K            = 100.00 * CCGO_EPSILON * 100 * 2
# Clementi Go angle force constant
const CCGO_ANGL_K            = 20.00 * CCGO_EPSILON * 2
# Clementi Go dihedral force constant
const CCGO_DIHE_K_1          = CCGO_EPSILON
const CCGO_DIHE_K_3          = CCGO_EPSILON * 0.5
# Clementi Go native contact eps
const CCGO_NATIVE_EPSILON    = CCGO_EPSILON

# ===============================
# Protein AICG2+ Model Parameters
# ===============================

# AICG2+ bond force constant
const AICG_BOND_K               = 110.40 * 100 * 2
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
const DNA3SPN_BOND_K_2    = 60.0 * 2 * JOU2CAL
# 3SPN.2C force constant for Gaussian dihedral
const DNA3SPN_DIH_G_K     = 7.0 * JOU2CAL
# 3SPN.2C sigma for Gaussian dihedral
const DNA3SPN_DIH_G_SIGMA = 0.3
# 3SPN.2C force constant for Gaussian dihedral
const DNA3SPN_DIH_P_K     = 2.0 * JOU2CAL

# ====================================
# RNA Structure-based Model Parameters
# ====================================

# RNA atomistic contact cutoff
const RNA_GO_ATOMIC_CUTOFF  = 5.5
# RNA stacking interaction dihedral cutoff
const RNA_STACK_DIH_CUTOFF  = 40.0
# RNA stacking interaction distance cutoff
const RNA_STACK_DIST_CUTOFF = 6.0
# RNA stacking interaction epsilon
const RNA_STACK_EPSILON     = 2.06
# RNA base pairing epsilon
const RNA_BPAIR_EPSILON_2HB = 2.94
const RNA_BPAIR_EPSILON_3HB = 5.37

RNA_BOND_K_LIST = Dict(
    "PS" => 26.5,
    "SR" => 40.3,
    "SY" => 62.9,
    "SP" => 84.1
)
RNA_ANGLE_K_LIST = Dict(
    "PSR" => 18.0,
    "PSY" => 22.8,
    "PSP" => 22.1,
    "SPS" => 47.8
)
RNA_DIHEDRAL_K_LIST = Dict(
    "PSPS" => 1.64,
    "SPSR" => 1.88,
    "SPSY" => 2.82,
    "SPSP" => 2.98
)
RNA_PAIR_EPSILON_OTHER = Dict(
    "SS" => 1.48,
    "BS" => 0.98,
    "SB" => 0.98,
    "BB" => 0.93
)

# =================
# PWMcos parameters
# =================
# PWMcos atomistic contact cutoff
const PWMCOS_ATOMIC_CUTOFF    = 4.0

# ======================
# Protein-RNA parameters
# ======================
# protein-RNA Go-term coefficient
const PRO_RNA_GO_EPSILON_B    = 0.62
const PRO_RNA_GO_EPSILON_S    = 0.74


# ====================
# GRO TOP File Options
# ====================

# "NREXCL" in "[moleculetype]"
const MOL_NR_EXCL             = 3
# "CGNR" in "[atoms]"
const AICG_ATOM_FUNC_NR       = 1
const DNA3SPN_ATOM_FUNC_NR    = 1
const RNA_ATOM_FUNC_NR        = 1
# "f" in "[bonds]"
const AICG_BOND_FUNC_TYPE     = 1
const CCGO_BOND_FUNC_TYPE     = 1
const DNA3SPN_BOND_FUNC2_TYPE = 1
const DNA3SPN_BOND_FUNC4_TYPE = 21
const RNA_BOND_FUNC_TYPE      = 1
# "f" in AICG-type "[angles]"
const AICG_ANG_G_FUNC_TYPE    = 21
# "f" in CCGO-type "[angles]"
const CCGO_ANG_FUNC_TYPE      = 1
# "f" in Flexible-type "[angles]"
const AICG_ANG_F_FUNC_TYPE    = 22
# "f" in DNA "[angles]"
const DNA3SPN_ANG_FUNC_TYPE   = 1
# "f" in RNA "[angles]"
const RNA_ANG_FUNC_TYPE       = 1
# "f" in AICG-type "[dihedral]"
const AICG_DIH_G_FUNC_TYPE    = 21
# "f" in CCGO-type "[dihedral]"
const CCGO_DIH_P_FUNC_TYPE    = 1
# "f" in Flexible-type "[dihedral]"
const AICG_DIH_F_FUNC_TYPE    = 22
# "f" in DNA Gaussian "[dihedral]"
const DNA3SPN_DIH_G_FUNC_TYPE = 21
# "f" in DNA Periodic "[dihedral]"
const DNA3SPN_DIH_P_FUNC_TYPE = 1
const DNA3SPN_DIH_P_FUNC_PERI = 1
# "f" in RNA Periodic "[dihedral]"
const RNA_DIH_FUNC_TYPE       = 1
# "f" in Go-contacts "[pairs]"
const AICG_CONTACT_FUNC_TYPE  = 2
# "f" in Go-contacts "[pairs]"
const CCGO_CONTACT_FUNC_TYPE  = 1
# "f" in RNA Go-contacts "[pairs]"
const RNA_CONTACT_FUNC_TYPE   = 2
# "f" in pro-RNA Go-contacts "[pairs]"
const RNP_CONTACT_FUNC_TYPE   = 2
# "f" in protein-DNA PWMcos "[pwmcos]"
const PWMCOS_FUNC_TYPE        = 1



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

function compute_vec_angle(vec1, vec2)
    n1 = norm(vec1)
    n2 = norm(vec2)
    return acos( dot(vec1, vec2) / n1 / n2) / pi * 180.0
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
    total_mass      = 0
    tmp_coor        = zeros(Float64, 3)
    for i in atom_indices
        a_mass      = ATOM_MASS_DICT[atom_names[i][1]]
        a_coor      = atom_coors[:, i]
        total_mass += a_mass
        tmp_coor   += a_coor * a_mass
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
        contact_count[AICG_ITYPE_LR_CT]  = 0
    end

    # control the number of salty bridge
    if contact_count[AICG_ITYPE_SS_SB]  >= 2
        contact_count[AICG_ITYPE_SS_QX] += contact_count[AICG_ITYPE_SS_SB] - 1
        contact_count[AICG_ITYPE_SS_SB]  = 1
    end

    return contact_count
end

# -----------------
# 3SPN.2C DNA model
# -----------------

function get_DNA3SPN_angle_param(angle_type, base_step)
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

# -------------------------
# RNA structure-based model
# -------------------------
function is_RNA_hydrogen_bond(atom_name_1, atom_name_2)
    special_atom_list = ['F', 'O', 'N']
    if atom_name_1 in special_atom_list && atom_name_2 in special_atom_list
        return true
    end
    return false
end

function compute_RNA_Go_contact(resid1, resid2, atom_names, atom_coors)
    hb_count = 0
    min_dist = 1e50
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

# ------------------------
# protein-DNA interactions
# ------------------------

function is_PWMcos_contact(resid1, resid2, atom_names, atom_coors)
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
            if dist_12 < PWMCOS_ATOMIC_CUTOFF
                return true
            end
        end
    end
    return false
end

# ------------------------
# protein-RNA interactions
# ------------------------

function is_protein_RNA_go_contact(resid1, resid2, atom_names, atom_coors)
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


# ------------------
# Other file formats
# ------------------

function read_modified_pfm(pfm_filename)
    pfm = Dict()
    for line in eachline(pfm_filename)
        words = split(line)
        if length(words) < 1
            continue
        end
        w1 = words[1]
        if occursin(w1, "ACGT")
            local_list = []
            for p in words[2:end]
                push!( local_list, parse(Float64, p) )
            end
            pfm[w1] = local_list
        elseif in(w1, ["CHAIN_A", "CHAIN_B"])
            local_list = []
            for dna_id in words[2:end]
                push!( local_list, parse(Int, dna_id) )
            end
            pfm[w1] = local_list
        end
    end

    pfmat = [pfm["A"]  pfm["C"]  pfm["G"]  pfm["T"]]
    ppmat = pfmat ./ sum(pfmat, dims=2)
    pwmat0 = -log.(ppmat)
    pwmat = pwmat0 .- sum(pwmat0, dims=2) ./ 4

    return (pwmat, pfm["CHAIN_A"], pfm["CHAIN_B"])
end


# =============================
# Coarse-Graining Structures!!!
# =============================

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
# core function
function pdb_2_top(args)

    # -----------------
    # Parsing arguments
    # -----------------
    pdb_name                = args["pdb"]
    protein_charge_filename = args["respac"]
    scale_scheme            = args["aicg-scale"]
    gen_3spn_itp            = args["3spn-param"]
    gen_pwmcos_itp          = args["pwmcos"]
    pwmcos_gamma            = args["pwmcos-scale"]
    pwmcos_epsil            = args["pwmcos-shift"]
    pfm_filename            = args["pfm"]
    appendto_filename       = args["patch"]
    do_output_psf           = args["psf"]
    do_output_cgpdb         = args["cgpdb"]
    do_debug                = args["debug"]
    do_output_sequence      = args["show-sequence"]
    ff_protein_name         = args["force-field-protein"]
    ff_DNA_name             = args["force-field-DNA"]
    ff_RNA_name             = args["force-field-RNA"]
    ccgo_contact_scale      = args["CCGO-contact-scale"]

    # ========================
    # Step -1: set force field
    # ========================
    if haskey(FF_PRO_DICT, ff_protein_name)
        ff_pro = FF_PRO_DICT[ff_protein_name]
    else
        error("Wrong force field for protein.")
    end
    if haskey(FF_DNA_DICT, ff_DNA_name)
        ff_dna = FF_DNA_DICT[ff_DNA_name]
    else
        error("Wrong force field for protein.")
    end
    if haskey(FF_RNA_DICT, ff_RNA_name)
        ff_rna = FF_RNA_DICT[ff_RNA_name]
    else
        error("Wrong force field for protein.")
    end

    # ===============
    # Step 0: numbers
    # ===============

    aa_num_atom    = 0
    aa_num_residue = 0
    aa_num_chain   = 0

    num_chain_pro  = 0
    num_chain_DNA  = 0
    num_chain_RNA  = 0

    i_step         = 0

    # ================
    # Step 1: open PDB
    # ================
    i_step += 1
    println("============================================================")
    println("> Step $(i_step): open PDB file.")

    aa_pdb_lines = []

    for line in eachline(pdb_name)
        if startswith(line, "ATOM")
            push!(aa_pdb_lines, rpad(line, 80))
            aa_num_atom += 1
        elseif startswith(line, "TER") || startswith(line, "END")
            push!(aa_pdb_lines, rpad(line, 80))
        end
    end

    aa_atom_name  = fill("    ",       aa_num_atom)
    aa_coor       = zeros(Float64, (3, aa_num_atom))

    aa_residues   = []
    aa_chains     = []

    i_atom        = 0
    i_resid       = 0
    curr_resid    = NaN
    curr_chain    = NaN
    curr_rname    = "    "
    residue_name  = "    "
    chain_id      = '?'
    tmp_res_atoms = []
    tmp_chain_res = []
    
    for line in aa_pdb_lines
        if startswith(line, "TER") || startswith(line, "END")
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

    println("          > Number of atoms    : $(aa_num_atom)")
    println("          > Number of residues : $(aa_num_residue)")
    println("          > Number of chains   : $(aa_num_chain)")

    # ===============================
    # Step 2: find out molecule types
    # ===============================
    i_step += 1
    println("============================================================")
    println("> Step $(i_step): set molecular types for every chain.")

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
                errmsg = @sprintf("BUG: Inconsistent residue types in chain %d ID - %s residue - %d : %s ",
                                  i_chain,
                                  chain.id,
                                  i_res,
                                  res_name)
                error(errmsg)
            end
        end
        cg_chain_mol_types[i_chain] = mol_type
        n_res = length(chain.residues)
        if mol_type == MOL_DNA
            n_particles = 3 * n_res - 1
            num_chain_DNA += 1
        elseif mol_type == MOL_RNA
            n_particles = 3 * n_res - 1
            num_chain_RNA += 1
        elseif mol_type == MOL_PROTEIN
            n_particles = n_res
            num_chain_pro += 1
        else
            n_particles = 0
        end
        cg_chain_length[i_chain] = n_particles
        cg_num_particles += n_particles
        @printf("          > Chain %3d | %7s \n", i_chain, MOL_TYPE_LIST[ mol_type ])
    end

    println("------------------------------------------------------------")
    @printf("          In total: %5d protein chains,\n", num_chain_pro)
    @printf("                    %5d DNA strands,\n", num_chain_DNA)
    @printf("                    %5d RNA strands.\n", num_chain_RNA)

    # ===========================
    # Step 3: Assign CG particles
    # ===========================
    i_step += 1
    println("============================================================")
    println("> Step $(i_step): assign coarse-grained particles.")

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
    cg_bead_type   = fill("    ", cg_num_particles)
    cg_bead_charge = zeros(Float64, cg_num_particles)
    cg_bead_mass   = zeros(Float64, cg_num_particles)
    cg_bead_coor   = zeros(Float64, (3, cg_num_particles))
    cg_chain_id    = zeros(Int, cg_num_particles)

    # protein
    top_cg_pro_bonds         = []
    top_cg_pro_angles        = []
    top_cg_pro_dihedrals     = []
    top_cg_pro_aicg13        = []
    top_cg_pro_aicg14        = []
    top_cg_pro_aicg_contact  = []

    param_cg_pro_e_13        = []
    param_cg_pro_e_14        = []
    param_cg_pro_e_contact   = []

    # DNA
    top_cg_DNA_bonds         = []
    top_cg_DNA_angles        = []
    top_cg_DNA_dih_Gaussian  = []
    top_cg_DNA_dih_periodic  = []

    # RNA
    top_cg_RNA_bonds         = []
    top_cg_RNA_angles        = []
    top_cg_RNA_dihedrals     = []
    top_cg_RNA_base_stack    = []
    top_cg_RNA_base_pair     = []
    top_cg_RNA_other_contact = []

    # protein-DNA
    top_cg_pro_DNA_pwmcos    = []

    # protein-RNA
    top_cg_pro_RNA_contact   = []

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

    if num_chain_pro > 0
        i_step += 1
        println("============================================================")
        println("> Step $(i_step): processing proteins.")

        # --------------------------------
        # Step 4.1: find out C-alpha atoms
        # --------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).1: determine CA mass, charge, and coordinates.")

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
                        cg_bead_type[i_res]    = res_name
                        cg_bead_charge[i_res]  = RES_CHARGE_DICT[res_name]
                        cg_bead_mass[i_res]    = RES_MASS_DICT[res_name]
                        cg_bead_coor[:, i_res] = aa_coor[:, i_atom]
                        cg_chain_id[i_res]     = i_chain
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
        println(">      $(i_step).2: AICG2+ topology.")
        println(" - - - - - - - - - - - - - - - - - - - - - - - -")
        println(">      $(i_step).2.1: AICG2+ local interactions.")
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
                coor1    = cg_bead_coor[:, i_res]
                coor2    = cg_bead_coor[:, i_res + 1]
                coor3    = cg_bead_coor[:, i_res + 2]
                dist13   = compute_distance(coor1, coor3)
                angle123 = compute_angle(coor1, coor2, coor3)
                push!(top_cg_pro_angles, ( i_res, angle123 ))
                push!(top_cg_pro_aicg13, ( i_res, dist13 ))

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
                coor1 = cg_bead_coor[:, i_res]
                coor2 = cg_bead_coor[:, i_res + 1]
                coor3 = cg_bead_coor[:, i_res + 2]
                coor4 = cg_bead_coor[:, i_res + 3]
                dihed = compute_dihedral(coor1, coor2, coor3, coor4)
                push!(top_cg_pro_dihedrals, ( i_res, dihed ))
                push!(top_cg_pro_aicg14, ( i_res, dihed ))

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
        println(">      $(i_step).2.2: AICG2+ Go-type native contacts.")
        e_ground_contact = 0.0
        num_contact = 0

        # intra-molecular contacts
        @printf("%11s Calculating intra-molecular contacts... \n", " ")
        @printf("              ... chain   : %32s", " ")
        i_progress_count = 0
        for i_chain in 1:aa_num_chain

            chain = cg_chains[i_chain]

            if chain.moltype != MOL_PROTEIN
                continue
            end

            # -----------------
            # show progress bar
            # -----------------
            i_progress_count += 1
            print("\b"^32)
            progress_percent = trunc(Int, i_progress_count / num_chain_pro * 20)
            progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
            @printf(" [%20s] %2d / %2d ", progress_bar, i_progress_count, num_chain_pro)
            # ------------------

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
        print("\n              ... intra-molecular contacts: DONE! \n")

        # inter-molecular ( protein-protein ) contacts
        if num_chain_pro > 1
            @printf("%11s Calculating inter-molecular contacts... \n", " ")
            @printf("              ... progress: %32s", " ")
            for i_chain in 1 : aa_num_chain - 1
                chain1 = cg_chains[i_chain]
    
                # -----------------
                # show progress bar
                # -----------------
                print("\b"^32)
                progress_percent = trunc(Int, i_chain / ( aa_num_chain - 1 ) * 20)
                progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / ( aa_num_chain - 1 ) * 100)
                # ------------------
    
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
            print("\n              ... inter-molecular contacts: DONE! \n")
        end

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

    end



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

    if num_chain_DNA > 0
        i_step += 1
        println("============================================================")
        println("> Step $(i_step): processing DNA.")

        # ----------------------------------
        #        Step 5.1: determine P, S, B
        # ----------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).1: determine P, S, B mass, charge, and coordinates.")

        for i_chain in 1 : aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_DNA
                continue
            end

            for i_res in chain.first : chain.last
                res_name  = cg_residues[i_res].res_name
                bead_name = cg_residues[i_res].atm_name
                bead_type = bead_name == "DP" || bead_name == "DS" ? bead_name : res_name
                bead_coor = compute_center_of_mass(cg_residues[i_res].atoms, aa_atom_name, aa_coor)
                cg_resid_name[i_res]   = res_name
                cg_resid_index[i_res]  = cg_residues[i_res].res_idx
                cg_bead_name[i_res]    = bead_name
                cg_bead_type[i_res]    = bead_type
                cg_bead_charge[i_res]  = RES_CHARGE_DICT[bead_type]
                cg_bead_mass[i_res]    = RES_MASS_DICT[bead_type]
                cg_bead_coor[:, i_res] = bead_coor[:]
                cg_chain_id[i_res]     = i_chain
            end
        end

        println(">           ... DONE!")

        # ---------------------------------
        #        Step 5.2: 3SPN.2C topology
        # ---------------------------------
        if gen_3spn_itp
            println("------------------------------------------------------------")
            println(">      $(i_step).2: 3SPN.2C topology.")
            println(" - - - - - - - - - - - - - - - - - - - - - - - -")
            println(">      $(i_step).2.1: 3SPN.2C local interactions.")
            for i_chain in 1:aa_num_chain

                chain = cg_chains[i_chain]

                if chain.moltype != MOL_DNA
                    continue
                end

                @printf("%11s Calculating DNA strand %d ... \n", " ", i_chain)
                @printf("              ... progress: %32s", " ")

                for i_res in chain.first : chain.last

                    # -----------------
                    # show progress bar
                    # -----------------
                    print("\b"^32)
                    progress_percent = trunc(Int, ( i_res - chain.first ) / ( chain.last - chain.first ) * 20)
                    progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                    @printf(" [%20s] %5.1f %% ", progress_bar, progress_percent * 5)
                    # ------------------

                    if cg_bead_name[i_res] == "DS"
                        # bond S--B
                        coor_s = cg_bead_coor[:, i_res]
                        coor_b = cg_bead_coor[:, i_res + 1]
                        r_sb   = compute_distance(coor_s, coor_b)
                        push!(top_cg_DNA_bonds, ( i_res, i_res + 1, r_sb ))
                        if i_res + 3 < chain.last
                            # bond S--P+1
                            coor_p3 = cg_bead_coor[:, i_res + 2]
                            r_sp3   = compute_distance(coor_s, coor_p3)
                            push!(top_cg_DNA_bonds, ( i_res, i_res + 2, r_sp3 ))
                            # Angle S--P+1--S+1
                            resname5  = cg_resid_name[i_res][end]
                            resname3  = cg_resid_name[i_res + 3][end]
                            coor_s3   = cg_bead_coor[:, i_res + 3]
                            ang_sp3s3 = compute_angle(coor_s, coor_p3, coor_s3)
                            k         = get_DNA3SPN_angle_param("SPS", resname5 * resname3)
                            push!(top_cg_DNA_angles, ( i_res, i_res + 2, i_res + 3, ang_sp3s3, k * 2 ))
                            # Dihedral S--P+1--S+1--B+1
                            coor_b3     = cg_bead_coor[:, i_res + 4]
                            dih_sp3s3b3 = compute_dihedral(coor_s, coor_p3, coor_s3, coor_b3)
                            push!(top_cg_DNA_dih_periodic, ( i_res, i_res + 2, i_res + 3, i_res + 4, dih_sp3s3b3 -180.0))
                            # Dihedral S--P+1--S+1--P+2
                            if i_res + 6 < chain.last
                                coor_p33     = cg_bead_coor[:, i_res + 5]
                                dih_sp3s3p33 = compute_dihedral(coor_s, coor_p3, coor_s3, coor_p33)
                                push!(top_cg_DNA_dih_periodic, ( i_res, i_res + 2, i_res + 3, i_res + 5, dih_sp3s3p33 - 180.0))
                                push!(top_cg_DNA_dih_Gaussian, ( i_res, i_res + 2, i_res + 3, i_res + 5, dih_sp3s3p33 ))
                            end
                        end
                    elseif cg_bead_name[i_res] == "DP"
                        # bond P--S
                        coor_p = cg_bead_coor[:, i_res]
                        coor_s = cg_bead_coor[:, i_res + 1]
                        r_ps   = compute_distance(coor_p, coor_s)
                        push!(top_cg_DNA_bonds, ( i_res, i_res + 1, r_ps ))
                        # angle P--S--B
                        resname5 = cg_resid_name[i_res - 1][end]
                        resname3 = cg_resid_name[i_res + 2][end]
                        coor_b   = cg_bead_coor[:, i_res + 2]
                        ang_psb  = compute_angle(coor_p, coor_s, coor_b)
                        k        = get_DNA3SPN_angle_param("PSB", resname5 * resname3)
                        push!(top_cg_DNA_angles, ( i_res, i_res + 1, i_res + 2, ang_psb, k * 2 ))
                        if i_res + 4 < chain.last
                            # angle P--S--P+1
                            coor_p3  = cg_bead_coor[:, i_res + 3]
                            ang_psp3 = compute_angle(coor_p, coor_s, coor_p3)
                            k        = get_DNA3SPN_angle_param("PSP", "all")
                            push!(top_cg_DNA_angles, ( i_res, i_res + 1, i_res + 3, ang_psp3, k * 2 ))
                            # Dihedral P--S--P+1--S+1
                            coor_s3    = cg_bead_coor[:, i_res + 4]
                            dih_psp3s3 = compute_dihedral(coor_p, coor_s, coor_p3, coor_s3)
                            push!(top_cg_DNA_dih_periodic, ( i_res, i_res + 1, i_res + 3, i_res + 4, dih_psp3s3 - 180.0))
                            push!(top_cg_DNA_dih_Gaussian, ( i_res, i_res + 1, i_res + 3, i_res + 4, dih_psp3s3 ))
                        end
                    elseif cg_bead_name[i_res] == "DB"
                        if i_res + 2 < chain.last
                            # angle B--S--P+1
                            resname5 = cg_resid_name[i_res][end]
                            resname3 = cg_resid_name[i_res + 1][end]
                            coor_b   = cg_bead_coor[:, i_res]
                            coor_s   = cg_bead_coor[:, i_res - 1]
                            coor_p3  = cg_bead_coor[:, i_res + 1]
                            ang_bsp3 = compute_angle(coor_b, coor_s, coor_p3)
                            k        = get_DNA3SPN_angle_param("BSP", resname5 * resname3)
                            push!(top_cg_DNA_angles, ( i_res, i_res - 1, i_res + 1, ang_bsp3, k * 2 ))
                            # Dihedral B--S--P+1--S+1
                            coor_s3    = cg_bead_coor[:, i_res + 2]
                            dih_bsp3s3 = compute_dihedral(coor_b, coor_s, coor_p3, coor_s3)
                            push!(top_cg_DNA_dih_periodic, ( i_res, i_res - 1, i_res + 1, i_res + 2, dih_bsp3s3 - 180.0))
                        end
                    else
                        errmsg = @sprintf("BUG: Wrong DNA particle type in chain %d, residue %d : %s ",
                                          i_chain,
                                          i_res,
                                          res_name)
                        error(errmsg)
                    end
                end
                print(" \n")
            end
            println(">           ... Bond, Angle, Dihedral: DONE!")
        end
    end


    # =========================
    # RNA structure based model
    # =========================
    #     ____  _   _    _    
    #    |  _ \| \ | |  / \   
    #    | |_) |  \| | / _ \  
    #    |  _ <| |\  |/ ___ \ 
    #    |_| \_\_| \_/_/   \_\
    # 
    # =========================

    if num_chain_RNA > 0
        i_step += 1
        println("============================================================")
        println("> Step $(i_step): processing RNA.")

        # ----------------------------------
        #         determine P, S, B
        # ----------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).1: determine P, S, B mass, charge, and coordinates.")

        for i_chain in 1 : aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_RNA
                continue
            end

            for i_res in chain.first : chain.last
                res_name  = cg_residues[i_res].res_name
                bead_name = cg_residues[i_res].atm_name
                bead_type = bead_name == "RP" || bead_name == "RS" ? bead_name : res_name
                # bead_coor = compute_center_of_mass(cg_residues[i_res].atoms, aa_atom_name, aa_coor)
                cg_resid_name[i_res]   = res_name
                cg_resid_index[i_res]  = cg_residues[i_res].res_idx
                cg_bead_name[i_res]    = bead_name
                cg_bead_type[i_res]    = bead_type
                cg_bead_charge[i_res]  = RES_CHARGE_DICT[bead_type]
                cg_bead_mass[i_res]    = RES_MASS_DICT[bead_type]
                cg_chain_id[i_res]     = i_chain
                if bead_name == "RP"
                    for i_atom in cg_residues[i_res].atoms
                        if aa_atom_name[i_atom][1] == 'P'
                            bead_coor = aa_coor[:, i_atom]
                        end
                    end
                elseif bead_name == "RS"
                    total_mass      = 0
                    tmp_coor        = zeros(Float64, 3)
                    for i_atom in cg_residues[i_res].atoms
                        a_name = aa_atom_name[i_atom]
                        if in(a_name, ["C1'", "C2'", "C3'", "C4'", "O4'"] )
                            a_mass      = ATOM_MASS_DICT[a_name[1]]
                            a_coor      = aa_coor[:, i_atom]
                            total_mass += a_mass
                            tmp_coor   += a_coor * a_mass
                        end
                    end
                    bead_coor = tmp_coor / total_mass
                elseif bead_name == "RB"
                    if res_name[end] == 'A' || res_name[end] == 'G'
                        for i_atom in cg_residues[i_res].atoms
                            if aa_atom_name[i_atom] == "N1"
                                bead_coor = aa_coor[:, i_atom]
                            end
                        end
                    else
                        for i_atom in cg_residues[i_res].atoms
                            if aa_atom_name[i_atom] == "N3"
                                bead_coor = aa_coor[:, i_atom]
                            end
                        end
                    end
                end
                cg_bead_coor[:, i_res] = bead_coor[:]
            end
        end

        println(">           ... DONE!")

        # -------------------------
        # Step 6.2: RNA topology
        # -------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).2: RNA topology.")
        println(" - - - - - - - - - - - - - - - - - - - - - - - -")
        println(">      $(i_step).2.1: RNA local interactions.")

        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_RNA
                continue
            end

            @printf("%11s Calculating bonded terms... \n", " ")
            for i_res in chain.first : chain.last
                if cg_bead_name[i_res] == "RS"
                    # bond S--B
                    coor_s    = cg_bead_coor[:, i_res]
                    coor_b    = cg_bead_coor[:, i_res + 1]
                    r_sb      = compute_distance(coor_s, coor_b)
                    base_type = cg_resid_name[i_res] in ["RA", "RG"] ? "R" : "Y"
                    bond_type = "S" * base_type
                    k         = RNA_BOND_K_LIST[bond_type]
                    push!(top_cg_RNA_bonds, (i_res, i_res + 1, r_sb , k * 2 * 100.0))
                    # bond S--P+1
                    if i_res + 2 < chain.last
                        coor_p3 = cg_bead_coor[:, i_res + 2]
                        r_sp3   = compute_distance(coor_s, coor_p3)
                        k       = RNA_BOND_K_LIST["SP"]
                        push!(top_cg_RNA_bonds, (i_res, i_res + 2, r_sp3 , k * 2 * 100.0))
                    end
                    if i_res + 4 <= chain.last
                        # Angle S--P+1--S+1
                        coor_s3   = cg_bead_coor[:, i_res + 3]
                        ang_sp3s3 = compute_angle(coor_s, coor_p3, coor_s3)
                        k         = RNA_ANGLE_K_LIST["SPS"]
                        push!(top_cg_RNA_angles, (i_res, i_res + 2, i_res + 3, ang_sp3s3, k * 2))
                        # Dihedral S--P+1--S+1--B+1
                        coor_b3     = cg_bead_coor[:, i_res + 4]
                        dih_sp3s3b3 = compute_dihedral(coor_s, coor_p3, coor_s3, coor_b3)
                        base_type   = cg_resid_name[i_res + 4] in ["RA", "RG"] ? "R" : "Y"
                        dihe_type   = "SPS" * base_type
                        k           = RNA_DIHEDRAL_K_LIST[dihe_type]
                        push!(top_cg_RNA_dihedrals, (i_res, i_res + 2, i_res + 3, i_res + 4, dih_sp3s3b3, k))
                    end
                    # Dihedral S--P+1--S+1--P+2
                    if i_res + 5 < chain.last
                        coor_p33     = cg_bead_coor[:, i_res + 5]
                        dih_sp3s3p33 = compute_dihedral(coor_s, coor_p3, coor_s3, coor_p33)
                        k            = RNA_DIHEDRAL_K_LIST["SPSP"]
                        push!(top_cg_RNA_dihedrals, (i_res, i_res + 2, i_res + 3, i_res + 5, dih_sp3s3p33, k))
                    end
                elseif cg_bead_name[i_res] == "RP"
                    # bond P--S
                    coor_p = cg_bead_coor[:, i_res]
                    coor_s = cg_bead_coor[:, i_res + 1]
                    r_ps   = compute_distance(coor_p, coor_s)
                    k      = RNA_BOND_K_LIST["PS"]
                    push!(top_cg_RNA_bonds, (i_res, i_res + 1, r_ps , k * 2 * 100.0))
                    # angle P--S--B
                    coor_b    = cg_bead_coor[:, i_res + 2]
                    ang_psb   = compute_angle(coor_p, coor_s, coor_b)
                    base_type = cg_resid_name[i_res + 2] in ["RA", "RG"] ? "R" : "Y"
                    angl_type = "PS" * base_type
                    k         = RNA_ANGLE_K_LIST[angl_type]
                    push!(top_cg_RNA_angles, (i_res, i_res + 1, i_res + 2, ang_psb, k * 2))
                    if i_res + 4 < chain.last
                        # angle P--S--P+1
                        coor_p3  = cg_bead_coor[:, i_res + 3]
                        ang_psp3 = compute_angle(coor_p, coor_s, coor_p3)
                        k        = RNA_ANGLE_K_LIST["PSP"]
                        push!(top_cg_RNA_angles, (i_res, i_res + 1, i_res + 3, ang_psp3, k * 2))
                        # Dihedral P--S--P+1--S+1
                        coor_s3    = cg_bead_coor[:, i_res + 4]
                        dih_psp3s3 = compute_dihedral(coor_p, coor_s, coor_p3, coor_s3)
                        k          = RNA_DIHEDRAL_K_LIST["PSPS"]
                        push!(top_cg_RNA_dihedrals, (i_res, i_res + 1, i_res + 3, i_res + 4, dih_psp3s3, k))
                    end
                elseif cg_bead_name[i_res] == "RB"
                    # do nothing...
                end
            end
        end

        # -----------------------
        # Go type native contacts
        # -----------------------
        println(" - - - - - - - - - - - - - - - - - - - - - - - -")
        println(">      $(i_step).2.2: RNA Go-type native contacts.")
        @printf("%11s Calculating intra-molecular contacts... \n", " ")
        @printf("              ... progress: %32s", " ")
        i_progress_count = 0
        for i_chain in 1:aa_num_chain
            chain = cg_chains[i_chain]

            if chain.moltype != MOL_RNA
                continue
            end

            # -----------------
            # show progress bar
            # -----------------
            i_progress_count += 1
            print("\b"^32)
            progress_percent = trunc(Int, i_progress_count / ( num_chain_RNA ) * 20)
            progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
            @printf(" [%20s] %2d / %2d ", progress_bar, i_progress_count, num_chain_RNA)
            # ------------------

            for i_res in chain.first : chain.last - 3

                if cg_bead_name[i_res] == "RP"
                    continue
                end

                coor_i = cg_bead_coor[:, i_res]

                for j_res in i_res + 3 : chain.last

                    if cg_bead_name[j_res] == "RP"
                        continue
                    end

                    if cg_bead_name[i_res] == "RS" || cg_bead_name[j_res] == "RS"
                        if j_res < i_res + 6
                            continue
                        end
                    end

                    coor_j = cg_bead_coor[:, j_res]
                    native_dist = compute_distance(coor_i, coor_j)
                    adist, nhb  = compute_RNA_Go_contact(cg_residues[i_res],
                                                         cg_residues[j_res],
                                                         aa_atom_name,
                                                         aa_coor)

                    if adist > RNA_GO_ATOMIC_CUTOFF
                        continue
                    end
                    
                    if j_res == i_res + 3 && cg_bead_name[i_res] == "RB"
                        coor_i_sug = cg_bead_coor[:, i_res - 1]
                        coor_j_sug = cg_bead_coor[:, j_res - 1]
                        st_dih = compute_dihedral(coor_i, coor_i_sug, coor_j_sug, coor_j)
                        if abs( st_dih ) < RNA_STACK_DIH_CUTOFF && adist < RNA_STACK_DIST_CUTOFF
                            push!(top_cg_RNA_base_stack, (i_res, j_res, native_dist, RNA_STACK_EPSILON))
                        else
                            push!(top_cg_RNA_other_contact, (i_res, j_res, native_dist, RNA_PAIR_EPSILON_OTHER["BB"]))
                        end
                    elseif cg_bead_name[i_res] == "RB" && cg_bead_name[j_res] == "RB"
                        if nhb == 2
                            push!(top_cg_RNA_base_pair, (i_res, j_res, native_dist, RNA_BPAIR_EPSILON_2HB))
                        elseif nhb >= 3
                            push!(top_cg_RNA_base_pair, (i_res, j_res, native_dist, RNA_BPAIR_EPSILON_3HB))
                        else
                            push!(top_cg_RNA_other_contact, (i_res, j_res, native_dist, RNA_PAIR_EPSILON_OTHER["BB"]))
                        end
                    else
                        contact_type = cg_bead_name[i_res][end] * cg_bead_name[j_res][end]
                        push!(top_cg_RNA_other_contact, (i_res, j_res, native_dist, RNA_PAIR_EPSILON_OTHER[contact_type]))
                    end
                end
            end
        end
        print("\n              ... intra-molecular contacts: DONE! \n")
 
        if num_chain_RNA > 1
            @printf("%11s Calculating inter-molecular contacts... \n", " ")
            @printf("              ... progress: %32s", " ")
            for i_chain in 1:aa_num_chain
    
                chain_1 = cg_chains[i_chain]
    
                if chain_1.moltype != MOL_RNA
                    continue
                end
    
                # -----------------
                # show progress bar
                # -----------------
                print("\b"^32)
                progress_percent = trunc(Int, i_chain / ( aa_num_chain - 1 ) * 20)
                progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
                @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / ( aa_num_chain - 1 ) * 100)
                # ------------------
    
                for i_res in chain_1.first : chain_1.last
                    if cg_bead_name[i_res] == "RP"
                        continue
                    end
                    coor_i = cg_bead_coor[:, i_res]
                    for j_chain in i_chain + 1 : aa_num_chain
                        chain_2 = cg_chains[j_chain]
                        if chain_2.moltype != MOL_RNA
                            continue
                        end
                        for j_res in chain_2.first : chain_2.last
                            if cg_bead_name[j_res] == "RP"
                                continue
                            end
                            coor_j = cg_bead_coor[:, j_res]
                            native_dist = compute_distance(coor_i, coor_j)
                            adist, nhb  = compute_RNA_Go_contact(cg_residues[i_res],
                                                                 cg_residues[j_res],
                                                                 aa_atom_name,
                                                                 aa_coor)
                            if adist > RNA_GO_ATOMIC_CUTOFF
                                continue
                            end
                            if cg_bead_name[i_res] == "RB" && cg_bead_name[j_res] == "RB"
                                if nhb == 2
                                    push!(top_cg_RNA_base_pair, (i_res, j_res, native_dist, RNA_BPAIR_EPSILON_2HB))
                                elseif nhb >= 3
                                    push!(top_cg_RNA_base_pair, (i_res, j_res, native_dist, RNA_BPAIR_EPSILON_3HB))
                                else
                                    push!(top_cg_RNA_other_contact, (i_res, j_res, native_dist, RNA_PAIR_EPSILON_OTHER["BB"]))
                                end
                            else
                                contact_type = cg_bead_name[i_res][end] * cg_bead_name[j_res][end]
                                push!(top_cg_RNA_other_contact, (i_res, j_res, native_dist, RNA_PAIR_EPSILON_OTHER[contact_type]))
                            end
                        end
                    end
                end
            end
            print("\n              ... inter-molecular contacts: DONE! \n")
        end
 
        println("------------------------------------------------------------")
        @printf("          > Total number of RNA contacts:     %12d  \n",
                length(top_cg_RNA_base_stack) + length(top_cg_RNA_base_pair) + length(top_cg_RNA_other_contact) )

    end


    # ===========================================================
    # Protein-RNA structure-based interactions: Go-like potential
    # ===========================================================
    #                  _       _             ____  _   _    _    
    #  _ __  _ __ ___ | |_ ___(_)_ __       |  _ \| \ | |  / \   
    # | '_ \| '__/ _ \| __/ _ \ | '_ \ _____| |_) |  \| | / _ \  
    # | |_) | | | (_) | ||  __/ | | | |_____|  _ <| |\  |/ ___ \ 
    # | .__/|_|  \___/ \__\___|_|_| |_|     |_| \_\_| \_/_/   \_\
    # |_|                                                        
    # 
    # ============================================================

    if num_chain_RNA > 0 && num_chain_pro > 0
        i_step += 1
        println("============================================================")
        println("> Step $(i_step): Generating protein-RNA native contacts.")

        @printf("%11s Calculating protein-RNA contacts... \n", " ")
        @printf("              ... progress: %32s", " ")
        for i_chain in 1:aa_num_chain
            chain_pro = cg_chains[i_chain]
            if chain_pro.moltype != MOL_PROTEIN
                continue
            end
            # -----------------
            # show progress bar
            # -----------------
            print("\b"^32)
            progress_percent = trunc(Int, i_chain / aa_num_chain * 20)
            progress_bar = "|" ^ progress_percent * " " ^ (20 - progress_percent)
            @printf(" [%20s] %5.1f %% ", progress_bar, i_chain / aa_num_chain * 100)
            # ------------------
    
            for i_res in chain_pro.first : chain_pro.last
                coor_i = cg_bead_coor[:, i_res]

                for j_chain in 1 : aa_num_chain
                    chain_RNA = cg_chains[j_chain]
                    if chain_RNA.moltype != MOL_RNA
                        continue
                    end
                    for j_res in chain_RNA.first : chain_RNA.last
                        if cg_bead_name[j_res] == "RP"
                            continue
                        end
                        if !is_protein_RNA_go_contact(cg_residues[i_res], cg_residues[j_res], aa_atom_name, aa_coor)
                            continue
                        end
                        coor_j = cg_bead_coor[:, j_res]
                        native_dist = compute_distance(coor_i, coor_j)
                        if cg_bead_name[j_res] == "RS"
                            push!(top_cg_pro_RNA_contact, (i_res, j_res, native_dist, PRO_RNA_GO_EPSILON_S))
                        elseif cg_bead_name[j_res] == "RB"
                            push!(top_cg_pro_RNA_contact, (i_res, j_res, native_dist, PRO_RNA_GO_EPSILON_B))
                        end
                    end
                end
            end
        end

        println(">           ... DONE!")
        println("------------------------------------------------------------")
        @printf("          > Total number of protein-RNA contacts: %8d  \n",
                length(top_cg_pro_RNA_contact) )
    end


    # ============================================================
    # PWMcos parameters: protein-DNA sequence-specific interaction
    # ============================================================
    #        ______        ____  __               
    #       |  _ \ \      / /  \/  | ___ ___  ___ 
    #       | |_) \ \ /\ / /| |\/| |/ __/ _ \/ __|
    #       |  __/ \ V  V / | |  | | (_| (_) \__ \
    #       |_|     \_/\_/  |_|  |_|\___\___/|___/
    # 
    # ============================================================
    if gen_pwmcos_itp
        pwmcos_native_contacts = []

        if num_chain_pro == 0
            error("Cannot generate PWMcos parameters without protein...")
        end
        if num_chain_DNA != 2
            error("Cannot generate PWMcos parameters from more or less than two DNA chains...")
        end

        i_step += 1
        println("============================================================")
        println("> Step $(i_step): Generating PWMcos parameters.")

        # ----------------------------------
        #        Step 7.1: determine P, S, B
        # ----------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).1: determine contacts between protein and DNA.")

        i_count_DNA = 0
        for i_chain in 1:aa_num_chain
            chain_pro = cg_chains[i_chain]

            if chain_pro.moltype != MOL_PROTEIN
                continue
            end

            for i_res in chain_pro.first : chain_pro.last
                i_res_N = i_res == chain_pro.first ? i_res : i_res - 1
                i_res_C = i_res == chain_pro.last  ? i_res : i_res + 1

                coor_pro_i = cg_bead_coor[:, i_res]
                coor_pro_N = cg_bead_coor[:, i_res_N]
                coor_pro_C = cg_bead_coor[:, i_res_C]

                for j_chain in 1:aa_num_chain
                    chain_DNA = cg_chains[j_chain]
                    
                    if chain_DNA.moltype != MOL_DNA
                        continue
                    end

                    for j_res in chain_DNA.first + 3 : chain_DNA.last - 3
                        if cg_bead_name[j_res] != "DB"
                            continue
                        end
                        if !is_PWMcos_contact(cg_residues[i_res], cg_residues[j_res], aa_atom_name, aa_coor)
                            continue
                        end

                        j_res_5, j_res_3 = j_res - 3, j_res + 3
                        coor_dna_j       = cg_bead_coor[:, j_res]
                        coor_dna_5       = cg_bead_coor[:, j_res_5]
                        coor_dna_3       = cg_bead_coor[:, j_res_3]
                        coor_dna_S       = cg_bead_coor[:, j_res - 1]

                        vec0   = coor_pro_i - coor_dna_j
                        vec1   = coor_dna_S - coor_dna_j
                        vec2   = coor_dna_3 - coor_dna_5
                        vec3   = coor_pro_N - coor_pro_C
                        r0     = norm(vec0)
                        theta1 = compute_vec_angle(vec0, vec1)
                        theta2 = compute_vec_angle(vec0, vec2)
                        theta3 = compute_vec_angle(vec0, vec3)

                        push!(pwmcos_native_contacts, (i_res - chain_pro.first + 1,
                                                       cg_resid_index[j_res],
                                                       r0,
                                                       theta1,
                                                       theta2,
                                                       theta3))

                        if do_debug
                            println("PWMcos | pro ===> ", i_res - chain_pro.first + 1,
                                    " DNA ===> ", j_res, " : ", cg_resid_index[j_res],
                                    " r0 = ", r0,
                                    " theta1 = ", theta1,
                                    " theta2 = ", theta2,
                                    " theta3 = ", theta3)
                        end
                    end
                end
            end
        end

        # ------------------------------------------------
        #        Step 7.2: Read in PFM and convert to PWM
        # ------------------------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).2: read in position frequency matrix (PFM).")

        pwmcos_pwm, pwmcos_chain_a, pwmcos_chain_b = read_modified_pfm(pfm_filename)
        num_pwmcos_terms = length(pwmcos_chain_a)


        # ------------------------------------------------
        #        Step 7.2: Read in PFM and convert to PWM
        # ------------------------------------------------
        println("------------------------------------------------------------")
        println(">      $(i_step).3: decomposing PWM.")

        ip_count = zeros(Float64, num_pwmcos_terms)

        contact_to_pwm = []
        for nat_contact in pwmcos_native_contacts
            i_dna = nat_contact[2] # cg_resid_index[dna]
            cnt_pwm_idx_a = indexin(i_dna, pwmcos_chain_a)[1]
            cnt_pwm_idx_b = indexin(i_dna, pwmcos_chain_b)[1]
            if cnt_pwm_idx_a isa Int
                push!( contact_to_pwm, (cnt_pwm_idx_a, 1) )
                ip_count[cnt_pwm_idx_a] += 1
            elseif cnt_pwm_idx_b isa Int
                push!( contact_to_pwm, (cnt_pwm_idx_b, -1) )
                ip_count[cnt_pwm_idx_b] += 1
            else
                error("Index error in CHAIN_A or CHAIN_B!")
            end
        end
        pwm_decomposed = pwmcos_pwm ./ ip_count

        for (i_cnt, nat_cnt) in enumerate(pwmcos_native_contacts)
            pwm_i, pwm_order = contact_to_pwm[i_cnt][1], contact_to_pwm[i_cnt][2]
            if pwm_order == 1
                eA, eC, eG, eT = pwm_decomposed[pwm_i, 1:4]
            elseif pwm_order == -1
                eA, eC, eG, eT = pwm_decomposed[pwm_i, 4:-1:1]
            end
            push!(top_cg_pro_DNA_pwmcos,
                  (nat_cnt[1],
                   nat_cnt[3],
                   nat_cnt[4],
                   nat_cnt[5],
                   nat_cnt[6],
                   eA, eC, eG, eT))
        end

        if do_debug
            println(size( contact_to_pwm ))
            println(pwm_decomposed)
        end

        println(">           ... DONE!")
    end




    # =========================================================================
    #   ___  _   _ _____ ____  _   _ _____ 
    #  / _ \| | | |_   _|  _ \| | | |_   _|
    # | | | | | | | | | | |_) | | | | | |  
    # | |_| | |_| | | | |  __/| |_| | | |  
    #  \___/ \___/  |_| |_|    \___/  |_|  
    # 
    # =========================================================================

    if gen_pwmcos_itp
        do_output_top    = false
        do_output_itp    = false
        do_output_gro    = false
        do_output_pwmcos = true
    else
        do_output_top    = true
        do_output_itp    = true
        do_output_gro    = true
        do_output_pwmcos = false
    end

    if do_output_sequence
        do_output_top    = false
        do_output_itp    = false
        do_output_gro    = false
        do_output_pwmcos = false
    end


    i_step += 1
    println("============================================================")
    println("> Step $(i_step): output .itp and .gro files.")

    # -------------------------------------------------------------------
    #        top ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # -------------------------------------------------------------------
    if do_output_top
        # ---------------
        #        filename
        # ---------------
        top_name = pdb_name[1:end-4] * "_cg.top"
        top_file = open(top_name, "w")

        itp_name = pdb_name[1:end-4] * "_cg.itp"
        itp_system_name = pdb_name[1:end-4]

        print(top_file, "; atom types for coarse-grained models\n")
        print(top_file, "#include \"./lib/atom_types.itp\" \n")
        if num_chain_pro > 0
            print(top_file, "; AICG2+ flexible local angle parameters \n")
            print(top_file, "#include \"./lib/flexible_local_angle.itp\" \n")
            print(top_file, "; AICG2+ flexible local dihedral parameters \n")
            print(top_file, "#include \"./lib/flexible_local_dihedral.itp\" \n")
        end
        print(top_file, "\n")

        print(top_file, "; Molecule topology \n")
        print(top_file, "#include \"./top/", itp_name, "\" \n\n")

        print(top_file, "[ system ] \n")
        print(top_file, itp_system_name, " \n\n")

        print(top_file, "[ molecules ] \n")
        print(top_file, itp_system_name, "  1 \n\n")

        print(top_file, "; [ cg_ele_mol_pairs ] \n")
        print(top_file, "; ON 1 - 2 : 3 - 4 \n")
        print(top_file, "; OFF 1 - 1 : 3 - 3 \n\n")

        print(top_file, "; [ pwmcos_mol_pairs ] \n")
        print(top_file, "; ON 1 - 2 : 3 - 4 \n")
        print(top_file, "; OFF 1 - 1 : 3 - 3 \n\n")
    end


    # -------------------------------------------------------------------
    #        itp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # -------------------------------------------------------------------
    if do_output_itp
        itp_mol_head     = "[ moleculetype ]\n"
        itp_mol_comm     = format(";{1:15s} {2:6s}\n", "name", "nrexcl")
        itp_mol_line     = "{1:<16} {2:>6d}\n"

        itp_atm_head     = "[ atoms ]\n"
        itp_atm_comm     = format(";{:>9}{:>5}{:>10}{:>5}{:>5}{:>5} {:>8} {:>8}\n", "nr", "type", "resnr", "res", "atom", "cg", "charge", "mass")
        itp_atm_line     = "{:>10d}{:>5}{:>10d}{:>5}{:>5}{:>5d} {:>8.3f} {:>8.3f}\n"

        itp_bnd_head     = "[ bonds ]\n"
        itp_bnd_comm     = format(";{:>9}{:>10}{:>5}{:>18}{:>18}\n", "i", "j", "f", "eq", "coef")
        itp_bnd_line     = "{:>10d}{:>10d}{:>5d}{:>18.4E}{:>18.4E}\n"

        itp_13_head      = "[ angles ] ; AICG2+ 1-3 interaction\n"
        itp_13_comm      = format(";{:>9}{:>10}{:>10}{:>5}{:>15}{:>15}{:>15}\n", "i", "j", "k", "f", "eq", "coef", "w")
        itp_13_line      = "{:>10d}{:>10d}{:>10d}{:>5d}{:>15.4E}{:>15.4E}{:>15.4E}\n"

        itp_ang_f_head   = "[ angles ] ; AICG2+ flexible local interaction\n"
        itp_ang_f_comm   = format(";{:>9}{:>10}{:>10}{:>5}\n", "i", "j", "k", "f")
        itp_ang_f_line   = "{:>10d}{:>10d}{:>10d}{:>5d}\n"

        itp_ang_head     = "[ angles ] ; cannonical angle \n"
        itp_ang_comm     = format(";{:>9}{:>10}{:>10}{:>5}{:>18}{:>18} \n", "i", "j", "k", "f", "eq", "coef")
        itp_ang_line     = "{:>10d}{:>10d}{:>10d}{:>5d}{:>18.4E}{:>18.4E}\n"

        itp_dih_P_head   = "[ dihedrals ] ; periodic dihedrals\n"
        itp_dih_P_comm   = format(";{:>9}{:>10}{:>10}{:>10}{:>5}{:>18}{:>18}{:>5}\n", "i", "j", "k", "l", "f", "eq", "coef", "n")
        itp_dih_P_line   = "{:>10d}{:>10d}{:>10d}{:>10d}{:>5d}{:>18.4E}{:>18.4E}{:>5d}\n"

        itp_dih_G_head   = "[ dihedrals ] ; Gaussian dihedrals\n"
        itp_dih_G_comm   = format(";{:>9}{:>10}{:>10}{:>10}{:>5}{:>15}{:>15}{:>15}\n", "i", "j", "k", "l", "f", "eq", "coef", "w")
        itp_dih_G_line   = "{:>10d}{:>10d}{:>10d}{:>10d}{:>5d}{:>15.4E}{:>15.4E}{:>15.4E}\n"

        itp_dih_F_head   = "[ dihedrals ] ; AICG2+ flexible local interation\n"
        itp_dih_F_comm   = format(";{:>9}{:>10}{:>10}{:>10}{:>5}\n", "i", "j", "k", "l", "f")
        itp_dih_F_line   = "{:>10d}{:>10d}{:>10d}{:>10d}{:>5d}\n"

        itp_contact_head = "[ pairs ] ; Go-type native contact\n"
        itp_contact_comm = format(";{:>9}{:>10}{:>10}{:>15}{:>15}\n", "i", "j", "f", "eq", "coef")
        itp_contact_line = "{:>10d}{:>10d}{:>10d}{:>15.4E}{:>15.4E}\n"

        itp_exc_head     = "[ exclusions ] ; Genesis exclusion list\n"
        itp_exc_comm     = format(";{:>9}{:>10}\n", "i", "j")
        itp_exc_line     = "{:>10d}{:>10d}\n"

        # ---------------
        #        filename
        # ---------------
        itp_name = pdb_name[1:end-4] * "_cg.itp"
        itp_file = open(itp_name, "w")

        # --------------------
        # Writing CG particles
        # --------------------
        # print molecule type information
        itp_system_name = pdb_name[1:end-4]
        print(itp_file, itp_mol_head)
        print(itp_file, itp_mol_comm)
        printfmt(itp_file, itp_mol_line, itp_system_name, MOL_NR_EXCL)
        print(itp_file,"\n")

        # -----------------------
        #               [ atoms ]
        # -----------------------

        print(itp_file, itp_atm_head)
        print(itp_file, itp_atm_comm)
        for i_bead in 1 : cg_num_particles
            printfmt(itp_file,
                     itp_atm_line,
                     i_bead,
                     cg_bead_type[i_bead],
                     cg_resid_index[i_bead],
                     cg_resid_name[i_bead],
                     cg_bead_name[i_bead],
                     AICG_ATOM_FUNC_NR,
                     cg_bead_charge[i_bead],
                     cg_bead_mass[i_bead])
        end
        print(itp_file,"\n")

        # ----------------
        #        [ bonds ]
        # ----------------

        if length(top_cg_pro_bonds) + length(top_cg_DNA_bonds) + length(top_cg_RNA_bonds) > 0
            print(itp_file, itp_bnd_head)
            print(itp_file, itp_bnd_comm)
            
            # AICG2+ bonds
            if ff_pro == FF_pro_AICG2p
                for i_bond in 1 : length(top_cg_pro_bonds)
                    printfmt(itp_file,
                             itp_bnd_line,
                             top_cg_pro_bonds[i_bond][1],
                             top_cg_pro_bonds[i_bond][1] + 1,
                             AICG_BOND_FUNC_TYPE,
                             top_cg_pro_bonds[i_bond][2] * 0.1,
                             AICG_BOND_K * CAL2JOU)
                end
            end

            # Clementi Go bonds
            if ff_pro == FF_pro_Clementi_Go
                for i_bond in 1 : length(top_cg_pro_bonds)
                    printfmt(itp_file,
                             itp_bnd_line,
                             top_cg_pro_bonds[i_bond][1],
                             top_cg_pro_bonds[i_bond][1] + 1,
                             CCGO_BOND_FUNC_TYPE,
                             top_cg_pro_bonds[i_bond][2] * 0.1,
                             CCGO_BOND_K * CAL2JOU)
                end
            end
            

            # 3SPN.2C bonds
            if ff_dna == FF_DNA_3SPN2C
                for i_bond in 1 : length(top_cg_DNA_bonds)
                    printfmt(itp_file,
                             itp_bnd_line,
                             top_cg_DNA_bonds[i_bond][1],
                             top_cg_DNA_bonds[i_bond][2],
                             DNA3SPN_BOND_FUNC4_TYPE,
                             top_cg_DNA_bonds[i_bond][3] * 0.1,
                             DNA3SPN_BOND_K_2 * CAL2JOU)
                end
            end

            # Structure-based RNA bonds
            if ff_rna == FF_RNA_Go
                for i_bond in 1 : length(top_cg_RNA_bonds)
                    printfmt(itp_file,
                             itp_bnd_line,
                             top_cg_RNA_bonds[i_bond][1],
                             top_cg_RNA_bonds[i_bond][2],
                             RNA_BOND_FUNC_TYPE,
                             top_cg_RNA_bonds[i_bond][3] * 0.1,
                             top_cg_RNA_bonds[i_bond][4] * CAL2JOU)
                end
            end

            print(itp_file, "\n")

        end

        # -----------------
        #        [ angles ]
        # -----------------

        # AICG2+ 1-3
        if ff_pro == FF_pro_AICG2p
            if length(top_cg_pro_aicg13) > 0
                print(itp_file, itp_13_head)
                print(itp_file, itp_13_comm)
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
                print(itp_file, "\n")
            end
        end

        # AICG2+ flexible
        if ff_pro == FF_pro_AICG2p
            if length(top_cg_pro_angles) > 0
                print(itp_file, itp_ang_f_head)
                print(itp_file, itp_ang_f_comm)
                for i_ang in 1 : length(top_cg_pro_angles)
                    printfmt(itp_file,
                             itp_ang_f_line,
                             top_cg_pro_angles[i_ang][1],
                             top_cg_pro_angles[i_ang][1] + 1,
                             top_cg_pro_angles[i_ang][1] + 2,
                             AICG_ANG_F_FUNC_TYPE)
                end
                print(itp_file, "\n")
            end
        end

        # Clementi Go angle
        if ff_pro == FF_pro_Clementi_Go
            if length(top_cg_pro_angles) > 0
                print(itp_file, itp_ang_head)
                print(itp_file, itp_ang_comm)
                for i_ang in 1 : length(top_cg_pro_angles)
                    printfmt(itp_file,
                             itp_ang_line,
                             top_cg_pro_angles[i_ang][1],
                             top_cg_pro_angles[i_ang][1] + 1,
                             top_cg_pro_angles[i_ang][1] + 2,
                             CCGO_ANG_FUNC_TYPE,
                             top_cg_pro_angles[i_ang][2],
                             CCGO_ANGL_K * CAL2JOU)
                end
                print(itp_file, "\n")
            end
        end

        # 3SPN.2C angles
        if ff_dna == FF_DNA_3SPN2C
            if length(top_cg_DNA_angles) > 0
                print(itp_file, itp_ang_head)
                print(itp_file, itp_ang_comm)
                for i_ang in 1 : length(top_cg_DNA_angles)
                    printfmt(itp_file,
                             itp_ang_line,
                             top_cg_DNA_angles[i_ang][1],
                             top_cg_DNA_angles[i_ang][2],
                             top_cg_DNA_angles[i_ang][3],
                             DNA3SPN_ANG_FUNC_TYPE,
                             top_cg_DNA_angles[i_ang][4],
                             top_cg_DNA_angles[i_ang][5] * CAL2JOU)
                end
                print(itp_file, "\n")
            end
        end

        # RNA structure-based angles
        if ff_rna == FF_RNA_Go
            if length(top_cg_RNA_angles) > 0
                print(itp_file, itp_ang_head)
                print(itp_file, itp_ang_comm)
                for i_ang in 1 : length(top_cg_RNA_angles)
                    printfmt(itp_file,
                             itp_ang_line,
                             top_cg_RNA_angles[i_ang][1],
                             top_cg_RNA_angles[i_ang][2],
                             top_cg_RNA_angles[i_ang][3],
                             RNA_ANG_FUNC_TYPE,
                             top_cg_RNA_angles[i_ang][4],
                             top_cg_RNA_angles[i_ang][5] * CAL2JOU)
                end
                print(itp_file, "\n")
            end
        end

        # --------------------
        #        [ dihedrals ]
        # --------------------

        if ff_pro == FF_pro_AICG2p
            # AICG2+ Gaussian dihedrals
            if length(top_cg_pro_aicg14) > 0
                print(itp_file, itp_dih_G_head)
                print(itp_file, itp_dih_G_comm)
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
                print(itp_file, "\n")
            end
    
            # AICG2+ flexible dihedrals
            if length(top_cg_pro_dihedrals) > 0
                print(itp_file, itp_dih_F_head)
                print(itp_file, itp_dih_F_comm)
                for i_dih in 1 : length(top_cg_pro_dihedrals)
                    printfmt(itp_file,
                             itp_dih_F_line,
                             top_cg_pro_dihedrals[i_dih][1],
                             top_cg_pro_dihedrals[i_dih][1] + 1,
                             top_cg_pro_dihedrals[i_dih][1] + 2,
                             top_cg_pro_dihedrals[i_dih][1] + 3,
                             AICG_DIH_F_FUNC_TYPE)
                end
                print(itp_file, "\n")
            end
        end

        # Clementi Go dihedral
        if ff_pro == FF_pro_Clementi_Go
            if length(top_cg_pro_dihedrals) > 0
                print(itp_file, itp_dih_P_head)
                print(itp_file, itp_dih_P_comm)
                for i_dih in 1 : length(top_cg_pro_dihedrals)
                    printfmt(itp_file,
                             itp_dih_P_line,
                             top_cg_pro_dihedrals[i_dih][1],
                             top_cg_pro_dihedrals[i_dih][1] + 1,
                             top_cg_pro_dihedrals[i_dih][1] + 2,
                             top_cg_pro_dihedrals[i_dih][1] + 3,
                             CCGO_DIH_P_FUNC_TYPE,
                             top_cg_pro_dihedrals[i_dih][2] - 180.0,
                             CCGO_DIHE_K_1 * CAL2JOU,
                             1)
                end
                for i_dih in 1 : length(top_cg_pro_dihedrals)
                    printfmt(itp_file,
                             itp_dih_P_line,
                             top_cg_pro_dihedrals[i_dih][1],
                             top_cg_pro_dihedrals[i_dih][1] + 1,
                             top_cg_pro_dihedrals[i_dih][1] + 2,
                             top_cg_pro_dihedrals[i_dih][1] + 3,
                             CCGO_DIH_P_FUNC_TYPE,
                             3.0 * top_cg_pro_dihedrals[i_dih][2] - 180.0,
                             CCGO_DIHE_K_3 * CAL2JOU,
                             3)
                end
                print(itp_file, "\n")
            end
        end

        if ff_dna == FF_DNA_3SPN2C
            # 3SPN.2C Gaussian dihedrals
            if length(top_cg_DNA_dih_Gaussian) > 0
                print(itp_file, itp_dih_G_head)
                print(itp_file, itp_dih_G_comm)
                for i_dih in 1 : length(top_cg_DNA_dih_Gaussian)
                    printfmt(itp_file,
                             itp_dih_G_line,
                             top_cg_DNA_dih_Gaussian[i_dih][1],
                             top_cg_DNA_dih_Gaussian[i_dih][2],
                             top_cg_DNA_dih_Gaussian[i_dih][3],
                             top_cg_DNA_dih_Gaussian[i_dih][4],
                             DNA3SPN_DIH_G_FUNC_TYPE,
                             top_cg_DNA_dih_Gaussian[i_dih][5],
                             DNA3SPN_DIH_G_K * CAL2JOU,
                             DNA3SPN_DIH_G_SIGMA)
                end
                print(itp_file, "\n")
            end
    
            # 3SPN.2C Periodic dihedrals
            if length(top_cg_DNA_dih_periodic) > 0
                print(itp_file, itp_dih_P_head)
                print(itp_file, itp_dih_P_comm)
                for i_dih in 1 : length(top_cg_DNA_dih_periodic)
                    printfmt(itp_file,
                             itp_dih_P_line,
                             top_cg_DNA_dih_periodic[i_dih][1],
                             top_cg_DNA_dih_periodic[i_dih][2],
                             top_cg_DNA_dih_periodic[i_dih][3],
                             top_cg_DNA_dih_periodic[i_dih][4],
                             DNA3SPN_DIH_P_FUNC_TYPE,
                             top_cg_DNA_dih_periodic[i_dih][5],
                             DNA3SPN_DIH_P_K * CAL2JOU,
                             DNA3SPN_DIH_P_FUNC_PERI)
                end
                print(itp_file, "\n")
            end
        end

        # RNA structure-based Periodic dihedrals
        if ff_rna == FF_RNA_Go
            if length(top_cg_RNA_dihedrals) > 0
                print(itp_file, itp_dih_P_head)
                print(itp_file, itp_dih_P_comm)
                for i_dih in 1 : length(top_cg_RNA_dihedrals)
                    printfmt(itp_file,
                             itp_dih_P_line,
                             top_cg_RNA_dihedrals[i_dih][1],
                             top_cg_RNA_dihedrals[i_dih][2],
                             top_cg_RNA_dihedrals[i_dih][3],
                             top_cg_RNA_dihedrals[i_dih][4],
                             RNA_DIH_FUNC_TYPE,
                             top_cg_RNA_dihedrals[i_dih][5] - 180.0,
                             top_cg_RNA_dihedrals[i_dih][6] * CAL2JOU,
                             1)
                end
                for i_dih in 1 : length(top_cg_RNA_dihedrals)
                    printfmt(itp_file,
                             itp_dih_P_line,
                             top_cg_RNA_dihedrals[i_dih][1],
                             top_cg_RNA_dihedrals[i_dih][2],
                             top_cg_RNA_dihedrals[i_dih][3],
                             top_cg_RNA_dihedrals[i_dih][4],
                             RNA_DIH_FUNC_TYPE,
                             3.0 * top_cg_RNA_dihedrals[i_dih][5] - 180.0,
                             top_cg_RNA_dihedrals[i_dih][6] / 2 * CAL2JOU,
                             3)
                end
                print(itp_file, "\n")
            end
        end

        # ----------------
        #        [ pairs ]
        # ----------------

        # print protein Go-type native contacts
        if ff_pro == FF_pro_AICG2p
            if length(top_cg_pro_aicg_contact) > 0
                print(itp_file, itp_contact_head)
                print(itp_file, itp_contact_comm)
                for i_c in 1 : length(top_cg_pro_aicg_contact)
                    printfmt(itp_file,
                             itp_contact_line,
                             top_cg_pro_aicg_contact[i_c][1],
                             top_cg_pro_aicg_contact[i_c][2],
                             AICG_CONTACT_FUNC_TYPE,
                             top_cg_pro_aicg_contact[i_c][3] * 0.1,
                             param_cg_pro_e_contact[i_c] * CAL2JOU)
                end
                print(itp_file, "\n")
            end
        end

        # Clementi Go native contacts
        if ff_pro == FF_pro_Clementi_Go
            if length(top_cg_pro_aicg_contact) > 0
                print(itp_file, itp_contact_head)
                print(itp_file, itp_contact_comm)
                for i_c in 1 : length(top_cg_pro_aicg_contact)
                    r = top_cg_pro_aicg_contact[i_c][3] * 0.1
                    v = 6.0 * CCGO_NATIVE_EPSILON * CAL2JOU * r^10 * ccgo_contact_scale
                    w = 5.0 * CCGO_NATIVE_EPSILON * CAL2JOU * r^12 * ccgo_contact_scale
                    printfmt(itp_file,
                             itp_contact_line,
                             top_cg_pro_aicg_contact[i_c][1],
                             top_cg_pro_aicg_contact[i_c][2],
                             CCGO_CONTACT_FUNC_TYPE,
                             v,
                             w)
                end
                print(itp_file, "\n")
            end
        end

        # print RNA Go-type native contacts
        if ff_rna == FF_RNA_Go
            if length(top_cg_RNA_base_stack) + length(top_cg_RNA_base_pair) + length(top_cg_RNA_other_contact) > 0
                print(itp_file, itp_contact_head)
                print(itp_file, itp_contact_comm)
                for i_c in 1 : length(top_cg_RNA_base_stack)
                    printfmt(itp_file,
                             itp_contact_line,
                             top_cg_RNA_base_stack[i_c][1],
                             top_cg_RNA_base_stack[i_c][2],
                             RNA_CONTACT_FUNC_TYPE,
                             top_cg_RNA_base_stack[i_c][3] * 0.1,
                             top_cg_RNA_base_stack[i_c][4] * CAL2JOU)
                end
                for i_c in 1 : length(top_cg_RNA_base_pair)
                    printfmt(itp_file,
                             itp_contact_line,
                             top_cg_RNA_base_pair[i_c][1],
                             top_cg_RNA_base_pair[i_c][2],
                             RNA_CONTACT_FUNC_TYPE,
                             top_cg_RNA_base_pair[i_c][3] * 0.1,
                             top_cg_RNA_base_pair[i_c][4] * CAL2JOU)
                end
                for i_c in 1 : length(top_cg_RNA_other_contact)
                    printfmt(itp_file,
                             itp_contact_line,
                             top_cg_RNA_other_contact[i_c][1],
                             top_cg_RNA_other_contact[i_c][2],
                             RNA_CONTACT_FUNC_TYPE,
                             top_cg_RNA_other_contact[i_c][3] * 0.1,
                             top_cg_RNA_other_contact[i_c][4] * CAL2JOU)
                end
                print(itp_file, "\n")
            end
        end


        # print protein-RNA native contacts
        if ff_pro == FF_pro_AICG2p && ff_rna == FF_RNA_Go
            if length(top_cg_pro_RNA_contact) > 0
                print(itp_file, itp_contact_head)
                print(itp_file, itp_contact_comm)
                for i_c in 1 : length(top_cg_pro_RNA_contact)
                    printfmt(itp_file,
                             itp_contact_line,
                             top_cg_pro_RNA_contact[i_c][1],
                             top_cg_pro_RNA_contact[i_c][2],
                             RNP_CONTACT_FUNC_TYPE,
                             top_cg_pro_RNA_contact[i_c][3] * 0.1,
                             top_cg_pro_RNA_contact[i_c][4] * CAL2JOU)
                end
                print(itp_file, "\n")
            end
        end


        # ---------------------
        #        [ exclusions ]
        # ---------------------

        # print Protein exclusion list
        if ff_pro == FF_pro_AICG2p || ff_pro == FF_pro_Clementi_Go
            if length(top_cg_pro_aicg_contact) > 0
                print(itp_file, itp_exc_head)
                print(itp_file, itp_exc_comm)
                for i_c in 1 : length(top_cg_pro_aicg_contact)
                    printfmt(itp_file,
                             itp_exc_line,
                             top_cg_pro_aicg_contact[i_c][1],
                             top_cg_pro_aicg_contact[i_c][2])
                end
                print(itp_file, "\n")
            end
        end

        # print RNA exclusion list
        if ff_rna == FF_RNA_Go
            if length(top_cg_RNA_base_stack) + length(top_cg_RNA_base_pair) + length(top_cg_RNA_other_contact) > 0
                print(itp_file, itp_exc_head)
                print(itp_file, itp_exc_comm)
                for i_c in 1 : length(top_cg_RNA_base_stack)
                    printfmt(itp_file,
                             itp_exc_line,
                             top_cg_RNA_base_stack[i_c][1],
                             top_cg_RNA_base_stack[i_c][2])
                end
                for i_c in 1 : length(top_cg_RNA_base_pair)
                    printfmt(itp_file,
                             itp_exc_line,
                             top_cg_RNA_base_pair[i_c][1],
                             top_cg_RNA_base_pair[i_c][2])
                end
                for i_c in 1 : length(top_cg_RNA_other_contact)
                    printfmt(itp_file,
                             itp_exc_line,
                             top_cg_RNA_other_contact[i_c][1],
                             top_cg_RNA_other_contact[i_c][2])
                end
                print(itp_file, "\n")
            end
        end


        # print protein-RNA exclusion contacts
        if ff_pro == FF_pro_AICG2p && ff_rna == FF_RNA_Go
            if length(top_cg_pro_RNA_contact) > 0
                print(itp_file, itp_exc_head)
                print(itp_file, itp_exc_comm)
                for i_c in 1 : length(top_cg_pro_RNA_contact)
                    printfmt(itp_file,
                             itp_exc_line,
                             top_cg_pro_RNA_contact[i_c][1],
                             top_cg_pro_RNA_contact[i_c][2])
                end
                print(itp_file, "\n")
            end
        end


        close(itp_file)
        println(">           ... .itp: DONE!")
       
    end

    # ------------------------------------------------------------
    # gro ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ------------------------------------------------------------
    if do_output_gro
        # HEAD: time in the unit of ps
        GRO_HEAD_STR  = "{}, t= {:>16.3f} \n"
        # ATOM NUM: free format int
        GRO_ATOM_NUM  = "{:>12d} \n"
        # XYZ: in the unit of nm!!!
        GRO_ATOM_LINE = "{:>5d}{:>5}{:>5}{:>5d} {:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f} \n"
        GRO_BOX_LINE  = "{:>15.4f}{:>15.4f}{:>15.4f} \n\n"

        gro_name = pdb_name[1:end-4] * "_cg.gro"
        gro_file = open(gro_name, "w")

        printfmt(gro_file, GRO_HEAD_STR, "CG model for GENESIS ", 0)
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
    end


    # ------------------------------------------------------------
    # PWMcos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ------------------------------------------------------------
    if do_output_pwmcos
        itp_pwmcos_head = "[ pwmcos ]\n"
        itp_pwmcos_comm     = format(";{:>5}{:>4}{:>9}{:>9}{:>9}{:>9}{:>12}{:>12}{:>12}{:>12}{:>8}{:>8}\n",
                                     "i", "f", "r0", "theta1", "theta2", "theta3",
                                     "ene_A", "ene_C", "ene_G", "ene_T",
                                     "gamma", "eps'")
        # itp_pwmcos_line = "{:>6d} {:>3d} {:>8.5f} {:>8.3f} {:>8.3f} {:>8.3f}{:>12.6f}{:>12.6f}{:>12.6f}{:>12.6f}{:>8.3f}{:>8.3f} \n"
        itp_pwmcos_line = "%6d %3d %8.5f %8.3f %8.3f %8.3f%12.6f%12.6f%12.6f%12.6f%8.3f%8.3f \n"

        if length( appendto_filename ) == 0
            itp_pwmcos_name = pdb_name[1:end-4] * "_cg_pwmcos.itp_patch"
            itp_pwmcos_file = open(itp_pwmcos_name, "w")
        else
            itp_pwmcos_name = appendto_filename
            itp_pwmcos_file = open(itp_pwmcos_name, "a")
        end

        print(itp_pwmcos_file, itp_pwmcos_head)
        print(itp_pwmcos_file, itp_pwmcos_comm)
        for itpterm in top_cg_pro_DNA_pwmcos
            @printf(itp_pwmcos_file,
                    "%6d %3d %8.5f %8.3f %8.3f %8.3f%12.6f%12.6f%12.6f%12.6f%8.3f%8.3f \n",
                    itpterm[1],
                    PWMCOS_FUNC_TYPE,
                    itpterm[2] * 0.1,
                    itpterm[3],
                    itpterm[4],
                    itpterm[5],
                    itpterm[6],
                    itpterm[7],
                    itpterm[8],
                    itpterm[9],
                    pwmcos_gamma,
                    pwmcos_epsil)
        end
        print(itp_pwmcos_file, "\n")

        close(itp_pwmcos_file)
        println(">           ... ", itp_pwmcos_name, " pwmcos.itp: DONE!")
    end

    # ------------------------------------------------------------
    # psf ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ------------------------------------------------------------
    if do_output_psf
        psf_head_str = "PSF CMAP \n\n"
        psf_title_str0 = "      3 !NTITLE \n"
        psf_title_str1 = "REMARKS PSF file created with Julia. \n"
        psf_title_str2 = "REMARKS System: {1}  \n"
        psf_title_str5 = "REMARKS ======================================== \n"
        psf_title_str6 = "       \n"
        psf_title_str = psf_title_str0 * psf_title_str1 * psf_title_str2 * psf_title_str5 * psf_title_str6
        psf_atom_title = " {:>6d} !NATOM \n"
        # PSF_ATOM_LINE = " {atom_ser} {seg_id} {res_ser} {res_name} {atom_name} {atom_type}  {charge}  {mass}   0"
        psf_atom_line = " {:>6d} {:>3} {:>5d} {:>3} {:>3} {:>5}  {:>10.6f}  {:>10.6f}          0 \n"
        chain_id_set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"

        psf_name = pdb_name[1:end-4] * "_cg.psf"
        psf_file = open(psf_name, "w")
        print(psf_file, psf_head_str)
        printfmt(psf_file, psf_title_str, pdb_name[1:end-4])
        printfmt(psf_file, psf_atom_title, cg_num_particles)
        for i_bead in 1 : cg_num_particles
            printfmt(psf_file,
                     psf_atom_line,
                     i_bead,
                     chain_id_set[cg_chain_id[i_bead]],
                     cg_resid_index[i_bead],
                     cg_resid_name[i_bead],
                     cg_bead_name[i_bead],
                     cg_bead_type[i_bead],
                     cg_bead_charge[i_bead],
                     cg_bead_mass[i_bead])
        end
        print(psf_file,"\n")

        close(psf_file)
        println(">           ... .psf: DONE!")
    end

    # ------------------------------------------------------------
    # cgpdb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ------------------------------------------------------------
    if do_output_cgpdb
        cg_pdb_name = pdb_name[1:end-4] * "_cg.pdb"
        cg_pdb_file = open(cg_pdb_name, "w")
        cg_pdb_atom_line = "ATOM  {:>5d} {:>4s}{:1}{:<4s}{:1}{:>4d}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}{:>10s}{:2s}{:2s} \n"
        chain_id_set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"
        tmp_chain_id = 0
        for i_bead in 1 : cg_num_particles
            if cg_chain_id[i_bead] > tmp_chain_id
                if tmp_chain_id > 0
                    print(cg_pdb_file, "TER\n")
                end
                tmp_chain_id = cg_chain_id[i_bead]
            end
            printfmt(cg_pdb_file,
                     cg_pdb_atom_line,
                     i_bead,
                     cg_bead_name[i_bead],
                     ' ',
                     cg_resid_name[i_bead],
                     chain_id_set[cg_chain_id[i_bead]],
                     cg_resid_index[i_bead],
                     ' ',
                     cg_bead_coor[1 , i_bead],
                     cg_bead_coor[2 , i_bead],
                     cg_bead_coor[3 , i_bead],
                     0.0,
                     0.0,
                     "",
                     "",
                     "")
        end
        print(cg_pdb_file,"TER\n")
        print(cg_pdb_file,"END\n")
        print(cg_pdb_file,"\n")

        close(cg_pdb_file)
        println(">           ... .pdb (CG) : DONE!")
    end

    # --------------------------------------------------------------
    # show sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # --------------------------------------------------------------
    if do_output_sequence
        chain_id_set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"
        cg_seq_name = pdb_name[1:end-4] * "_cg.fasta"
        cg_seq_file = open(cg_seq_name, "w")

        for i_chain in 1:aa_num_chain
            chain = aa_chains[i_chain]
            mol_type = cg_chain_mol_types[i_chain]
            printfmt(cg_seq_file,
                     "> Chain {1} : {2} \n",
                     chain_id_set[i_chain],
                     MOL_TYPE_LIST[mol_type])

            for i_res in chain.residues
                res_name = aa_residues[i_res].name
                print(cg_seq_file, RES_SHORTNAME_DICT[res_name])
            end

            print(cg_seq_file, "\n")
        end

        close(cg_seq_file)
        println(">           ... sequence output : DONE!")
    end

    println("------------------------------------------------------------")
    println("------------------------------------------------------------")
    println("[1;32m FINISH! [0m ")
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

        "--force-field-protein"
        help = "Force field for protein."
        arg_type = String
        default = "AICG2+"

        "--force-field-DNA"
        help = "Force field for DNA."
        arg_type = String
        default = "3SPN.2C"

        "--force-field-RNA"
        help = "Force field for RNA."
        arg_type = String
        default = "Go"

        "--CCGO-contact-scale"
        help = "Scaling native contact interaction coefficient."
        arg_type = Float64
        default = 1.0

        "--respac", "-c"
        help = "RESPAC protein charge distribution data."
        arg_type = String
        default = ""

        "--aicg-scale"
        help = "Scale AICG2+ local interactions: 0) average; 1) general (default)."
        arg_type = Int
        default = 1

        "--3spn-param"
        help = "Generate 3SPN.2C parameters from x3DNA generated PDB structure."
        action = :store_true

        "--pwmcos"
        help = "Generate parameters for protein-DNA sequence-specific interactions."
        action = :store_true

        "--pwmcos-scale"
        help = "Energy scaling factor for PWMcos."
        arg_type = Float64
        default = 1.0

        "--pwmcos-shift"
        help = "Energy shifting factor for PWMcos."
        arg_type = Float64
        default = 0.0

        "--psf"
        help = "Prepare PSF file."
        action = :store_true

        "--cgpdb"
        help = "Prepare CG PDB file."
        action = :store_true

        "--pfm", "-p"
        help = "Position frequency matrix file for protein-DNA sequence-specific interactions."
        arg_type = String
        default = ""

        "--patch"
        help = "Append (apply patch) to .itp file."
        arg_type = String
        default = ""

        "--show-sequence"
        help = "Show sequence of molecules in PDB."
        action = :store_true

        "--debug"
        help = "DEBUG."
        action = :store_true
    end

    return parse_args(s)
end

# ====
# Main
# ====

function main()

    args = parse_commandline()

    pdb_2_top(args)

end

main()
