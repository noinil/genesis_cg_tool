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
    "VAL" =>  99.14
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
    "VAL" =>  0.0
)

RES_NAME_LIST_PROTEIN = (
    "ALA", "ARG", "ASN", "ASP",
    "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    "CYM", "CYT")

# RES_NAME_LIST_DNA = ("DA", "DC", "DG", "DT", "ADE", "GUA", "URA", "CYT")
RES_NAME_LIST_DNA = ("DA", "DC", "DG", "DT")

# RES_NAME_LIST_RNA = ("RA", "RC", "RG", "RT", "ADE", "GUA", "THY", "CYT")
RES_NAME_LIST_RNA = ("RA", "RC", "RG", "RU")

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
const AICG_ANG_GAUSS_SIGMA      = 0.15 * 0.1  # nm
# AICG2+ sigma for Gaussian dihedral
const AICG_DIH_GAUSS_SIGMA      = 0.15        # Rad ??
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

AICG_PAIRWISE_ENERGY = zeros(Float64, (17, ))
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


# ====================
# GRO TOP File Options
# ====================

# "NREXCL" in "[moleculetype]"
const MOL_NR_EXCL            = 3
# "CGNR" in "[atoms]"
const AICG_ATOM_FUNC_NR      = 1
# "f" in "[bonds]"
const AICG_BOND_FUNC_TYPE    = 1
# "f" in AICG-type "[angles]"
const AICG_ANG_G_FUNC_TYPE   = 21
# "f" in Flexible-type "[angles]"
const AICG_ANG_F_FUNC_TYPE   = 22
# "f" in AICG-type "[dihedral]"
const AICG_DIH_G_FUNC_TYPE   = 21
# "f" in Flexible-type "[dihedral]"
const AICG_DIH_F_FUNC_TYPE   = 22
# "f" in Go-contacts "[pairs]"
const AICG_CONTACT_FUNC_TYPE = 2


###############################################################################
#                                  Functions                                  #
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
            if dist_12 < AICG_GO_ATOMIC_CUTOFF
                return True
            end
        end
    end
    return False
end
function count_aicg_atomic_contact(resid1, resid2, res_name_1, res_name_2, atom_names, atom_coors)
    contact_count                   = [0 for _ in 1:17]
    contact_count[AICG_ITYPE_OFFST] = 1
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
            if dist_12 < AICG_GO_ATOMIC_CUTOFF
                contact_count[AICG_ITYPE_LR_CT] += 1
            elseif dist_12 < AICG_ATOMIC_CUTOFF
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

function pdb_2_top(pdb_name)

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
    # aa_resid_index = zeros(Int,         aa_num_atom)
    # aa_resid_name  = fill("    ",       aa_num_atom)
    # aa_chain_id    = fill('?',          aa_num_atom)
    aa_coor        = zeros(Float64, (3, aa_num_atom))
    # aa_occupancy   = zeros(Float64,     aa_num_atom)
    # aa_tempfactor  = zeros(Float64,     aa_num_atom)
    # aa_segment_id  = fill("          ", aa_num_atom)
    # aa_element     = fill("  ",         aa_num_atom)
    # aa_charge      = zeros(Float64,     aa_num_atom)

    aa_residues = []
    aa_chains   = []

    i_atom     = 0
    i_resid    = 0
    curr_resid = NaN
    curr_chain = NaN
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
        # aa_resid_index[i_atom] = residue_serial
        # aa_resid_name[i_atom]  = residue_name
        # aa_chain_id[i_atom]    = chain_id
        aa_coor[1, i_atom]     = coor_x
        aa_coor[2, i_atom]     = coor_y
        aa_coor[3, i_atom]     = coor_z
        # aa_occupancy[i_atom]   = occupancy
        # aa_tempfactor[i_atom]  = tempfactor
        # aa_segment_id[i_atom]  = segment_id
        # aa_element[i_atom]     = element_name
        # aa_charge[i_atom]      = charge

        if residue_serial != curr_resid
            i_resid += 1
            push!(tmp_chain_res, i_resid)
            curr_resid = residue_serial
            if length(tmp_res_atoms) > 0
                push!(aa_residues, AAResidue(residue_name, tmp_res_atoms))
                tmp_res_atoms = []
            end
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
    
    aa_chain_mol_types = ones(Int, (aa_num_chain, ))
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
        aa_chain_mol_types[i_chain] = mol_type
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
        cg_num_particles += n_particles
        @printf("          > Chain %3d | %7s | # particles: %5d \n",
                i_chain, MOL_TYPE_LIST[ mol_type ], n_particles)
    end

    # =================================
    # Step 3: AICG2+ model for proteins
    # =================================
    println("============================================================")
    println("> Step 3: processing proteins.")

    # --------------------------------
    # Step 3.1: find out C-alpha atoms
    # --------------------------------


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
    end
    return parse_args(s)
end

# ====
# Main
# ====

function main()
    
    args = parse_commandline()

    pdb_2_top(args["pdb"])

end

main()

