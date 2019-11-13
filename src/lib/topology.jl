###############################################################################
#                                   Topology                                  #
###############################################################################


# ========================================
# Coarse Grained Model Topology Structures
# ========================================

struct CGTopBond
    i::Int
    j::Int
    r0::Float64
end

struct CGTopAngle
    i::Int
    j::Int
    k::Int
    a0::Float64
end

struct CGTopDihedral
    i::Int
    j::Int
    k::Int
    l::Int
    t0::Float64
end

struct CGTopContact
    i::Int
    j::Int
    r0::Float64
end

struct CGTopPWMcos
    i::Int
    r0::Float64
    t1::Float64
    t2::Float64
    t3::Float64
    eA::Float64
    eC::Float64
    eG::Float64
    eT::Float64
end

struct CGTopology
    # CG particles
    cg_resid_name::Vector{String}
    cg_resid_index::Vector{Int}
    cg_bead_name::Vector{String}
    cg_bead_type::Vector{String}
    cg_bead_charge::Vector{<:Real}
    cg_bead_mass::Vector{<:Real}
    cg_chain_id::Vector{Int}
    cg_seg_name::Vector{String}
    # protein
    top_cg_pro_bonds::Vector{CGTopBond}
    top_cg_pro_angles::Vector{CGTopAngle}
    top_cg_pro_dihedrals::Vector{CGTopDihedral}
    top_cg_pro_aicg13::Vector{CGTopAngle}
    top_cg_pro_aicg14::Vector{CGTopDihedral}
    top_cg_pro_go_contact::Vector{CGTopContact}
    param_cg_pro_e_13::Vector{Float64}
    param_cg_pro_e_14::Vector{Float64}
    param_cg_pro_e_contact::Vector{Float64}
    # DNA
    top_cg_DNA_bonds::Vector{CGTopBond}
    top_cg_DNA_angles::Vector{CGTopAngle}
    top_cg_DNA_dih_Gaussian::Vector{CGTopDihedral}
    top_cg_DNA_dih_periodic::Vector{CGTopDihedral}
    param_cg_DNA_k_angles::Vector{Float64}
    # RNA
    top_cg_RNA_bonds::Vector{CGTopBond}
    top_cg_RNA_angles::Vector{CGTopAngle}
    top_cg_RNA_dihedrals::Vector{CGTopDihedral}
    top_cg_RNA_base_stack::Vector{CGTopContact}
    top_cg_RNA_base_pair::Vector{CGTopContact}
    top_cg_RNA_other_contact::Vector{CGTopContact}
    param_cg_RNA_k_bonds::Vector{Float64}
    param_cg_RNA_k_angles::Vector{Float64}
    param_cg_RNA_k_dihedrals::Vector{Float64}
    param_cg_RNA_e_base_stack::Vector{Float64}
    param_cg_RNA_e_base_pair::Vector{Float64}
    param_cg_RNA_e_other_contact::Vector{Float64}
    # protein-DNA
    top_cg_pro_DNA_pwmcos::Vector{CGTopPWMcos}
    # protein-RNA
    top_cg_pro_RNA_contact::Vector{CGTopContact}
    param_cg_pro_RNA_e_contact::Vector{Float64}
end





# ===========================
# General Topology Structures
# ===========================

