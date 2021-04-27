###############################################################################
#                                  Parameters                                 #
###############################################################################

# ==================
# Physical Constants
# ==================
const CAL2JOU = 4.184
const JOU2CAL = 1.0 / CAL2JOU

# ============
# Force fields
# ============

struct ForceFieldCG
    ff_protein::Int
    ff_DNA::Int
    ff_RNA::Int
    ff_protein_DNA::Int
    ff_protein_RNA::Int
    ff_DNA_RNA::Int
end

# protein
const FF_pro_AICG2p      = 1
const FF_pro_Clementi_Go = 2
const FF_pro_KB_Go       = 3
# DNA
const FF_DNA_3SPN2C      = 1
# RNA
const FF_RNA_HT          = 1
# protein-DNA
const FF_PWMcos          = 1
const FF_pro_DNA_Go      = 2
const FF_PWMcos_ns       = 3
# protein-RNA
const FF_pro_RNA_Go      = 1
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
    "HT"       => FF_RNA_HT
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
    "HSD" => 137.14,
    "HSE" => 137.14,
    "HSP" => 138.14,
    "HID" => 137.14,
    "HIE" => 137.14,
    "HIP" => 138.14,
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
    "RT"  => 125.10,
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
    "HSD" =>  0.0,
    "HSE" =>  0.0,
    "HSP" =>  1.0,
    "HID" =>  0.0,
    "HIE" =>  0.0,
    "HIP" =>  1.0,
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
    "RT"  =>  0.0,
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
    "HSD" => "H",
    "HSE" => "H",
    "HSP" => "H",
    "HID" => "H",
    "HIE" => "H",
    "HIP" => "H",
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
    "RU"  => "U",
    "ADE" => "A",
    "CYT" => "C",
    "GUA" => "G",
    "URA" => "U",
    "THY" => "T",
    "A"   => "A",
    "C"   => "C",
    "G"   => "G",
    "U"   => "U",
    "T"   => "T"
)

RES_NAME_LIST_PROTEIN = (
    "ALA", "ARG", "ASN", "ASP",
    "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    "HSD", "HSE", "HSP",
    "HID", "HIE", "HIP",
    "CYM", "UNK")
RES_NAME_PROTEIN_DICT = Dict(
    "ALA" => "ALA",
    "ARG" => "ARG",
    "ASN" => "ASN",
    "ASP" => "ASP",
    "CYS" => "CYS",
    "CYM" => "CYS",
    "GLN" => "GLN",
    "GLU" => "GLU",
    "GLY" => "GLY",
    "HIS" => "HIS",
    "HSD" => "HIS",
    "HSE" => "HIS",
    "HSP" => "HIS",
    "HID" => "HIS",
    "HIE" => "HIS",
    "HIP" => "HIS",
    "ILE" => "ILE",
    "LEU" => "LEU",
    "LYS" => "LYS",
    "MET" => "MET",
    "PHE" => "PHE",
    "PRO" => "PRO",
    "SER" => "SER",
    "THR" => "THR",
    "TRP" => "TRP",
    "TYR" => "TYR",
    "VAL" => "VAL",
    "UNK" => "UNK"
)

RES_NAME_LIST_DNA = ("DA", "DC", "DG", "DT")
RES_NAME_DNA_DICT = Dict(
    "DA"  => "DA",
    "DC"  => "DC",
    "DG"  => "DG",
    "DT"  => "DT",
    "A"   => "DA",
    "C"   => "DC",
    "G"   => "DG",
    "T"   => "DT",
    "ADE" => "DA",
    "CYT" => "DC",
    "GUA" => "DG",
    "THY" => "DT"
)

RES_NAME_LIST_RNA = ("RA", "RC", "RG", "RU")
RES_NAME_RNA_DICT = Dict(
    "RA"  => "RA",
    "RC"  => "RC",
    "RG"  => "RG",
    "RT"  => "RT",
    "RU"  => "RU",
    "A"   => "RA",
    "C"   => "RC",
    "G"   => "RG",
    "T"   => "RT",
    "U"   => "RU",
    "ADE" => "RA",
    "CYT" => "RC",
    "GUA" => "RG",
    "URA" => "RU",
    "THY" => "RT"
)


# DNA CG residue atom names
ATOM_NAME_LIST_DP = ("P", "OP1", "OP2", "O5'", "O1P", "O2P")
ATOM_NAME_LIST_DS = ("C5'", "C4'", "C3'", "C2'", "C1'", "O4'")

# RNA CG residue atom names
ATOM_NAME_LIST_RP = ("P", "OP1", "OP2", "O1P", "O2P")
ATOM_NAME_LIST_RS = ("C5'", "C4'", "C3'", "C2'", "C1'", "O5'", "O4'", "O3'", "O2'")

RES_FASTA_LONGNAME_DICT_PRO = Dict(
    'A' => "ALA",
    'R' => "ARG",
    'N' => "ASN",
    'D' => "ASP",
    'C' => "CYS",
    'Q' => "GLN",
    'E' => "GLU",
    'G' => "GLY",
    'H' => "HIS",
    'I' => "ILE",
    'L' => "LEU",
    'K' => "LYS",
    'M' => "MET",
    'F' => "PHE",
    'P' => "PRO",
    'S' => "SER",
    'T' => "THR",
    'W' => "TRP",
    'Y' => "TYR",
    'V' => "VAL"
)


# ==============
# Molecule Types
# ==============

const MOL_DNA     = 1
const MOL_RNA     = 2
const MOL_PROTEIN = 3
const MOL_OTHER   = 4
MOL_TYPE_LIST = ("DNA", "RNA", "protein", "other", "unknown")

# ===========================
# General thresholds, cutoffs
# ===========================

const CG_MOL_CONTACT_CUTOFF = 20.0

const DIHEDRAL_SAFE_CUTOFF = 150.0
const DIHEDRAL_GAUS_MOD_TYPE = Dict(
    0 => 21,                    # use-dafe-dihedral = 0
    1 => 41,                    # use-safe-dihedral = 1; cos^2(kθ) type
    3 => 43                     # use-safe-dihedral = 3; sin^3(kθ) type
)
const DIHEDRAL_PERI_MOD_TYPE = Dict(
    0 => 1,                     # use-dafe-dihedral = 0
    1 => 32,                    # use-safe-dihedral = 1; cos^2(kθ) type
    2 => 31,                    # use-safe-dihedral = 2; remove dangerous dih
    3 => 33                     # use-safe-dihedral = 3; sin^3(θ) type
)
const DIHEDRAL_TABU_MOD_TYPE = Dict(
    0 => 22,                    # use-dafe-dihedral = 0
    1 => 52                     # use-safe-dihedral = 1; cos^2(kθ) type
)

###############################################################################
#                         Molecule specific parameters                        #
###############################################################################

# ====================================
# Protein Clementi Go Model Parameters
# ====================================

# Clementi Go energy unit: epsilon
const CCGO_EPSILON           = 1.0
# Clementi Go bond force constant
const CCGO_BOND_K            = 100.00 * CCGO_EPSILON
# Clementi Go angle force constant
const CCGO_ANGL_K            = 20.00 * CCGO_EPSILON
# Clementi Go dihedral force constant
const CCGO_DIHE_K_1          = CCGO_EPSILON
const CCGO_DIHE_K_3          = CCGO_EPSILON * 0.5
# Clementi Go native contact eps
const CCGO_NATIVE_EPSILON    = CCGO_EPSILON

# ===============================
# Protein AICG2+ Model Parameters
# ===============================

# AICG2+ bond force constant
const AICG_BOND_K               = 110.40
# AICG2+ sigma for Gaussian angle
const AICG_13_SIGMA             = 0.15  # A
# AICG2+ sigma for Gaussian dihedral
const AICG_14_SIGMA             = 0.15  # Rad ??
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
const DNA3SPN_BOND_K_2    = 0.6 * JOU2CAL
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
    "BB" => 0.93,
    # undetermined...
    "PP" => 1.00,
    "PS" => 1.00,
    "SP" => 1.00,
    "PB" => 1.00,
    "BP" => 1.00
)

# =================
# PWMcos parameters
# =================
# PWMcos atomistic contact cutoff
const PWMCOS_ATOMIC_CUTOFF     = 4.0
const pro_DNA_GO_ATOMIC_CUTOFF = 6.5

# ======================
# Protein-RNA parameters
# ======================
# protein-RNA Go-term coefficient
const PRO_RNA_GO_EPSILON_B    = 0.62
const PRO_RNA_GO_EPSILON_S    = 0.74
const PRO_RNA_GO_EPSILON_P    = 0.50


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
const CCGO_CONTACT_FUNC_TYPE  = 2
# "f" in RNA Go-contacts "[pairs]"
const RNA_CONTACT_FUNC_TYPE   = 2
# "f" in pro-RNA Go-contacts "[pairs]"
const RNP_CONTACT_FUNC_TYPE   = 2
# "f" in protein-DNA PWMcos "[pwmcos]"
const PWMCOS_FUNC_TYPE        = 1
const PWMCOS_NS_FUNC_TYPE     = 2


