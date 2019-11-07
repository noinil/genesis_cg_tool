###############################################################################
#                                   Topology                                  #
###############################################################################

struct TopBond
    i::Int
    j::Int
    r0::Float64
end

struct TopAngle
    i::Int
    j::Int
    k::Int
    a0::Float64
end

struct TopDihedral
    i::Int
    j::Int
    k::Int
    l::Int
    t0::Float64
end

struct TopContact
    i::Int
    j::Int
    r0::Float64
end

struct TopPWMcos
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
    top_cg_pro_bonds::Vector{TopBond}
    top_cg_pro_angles::Vector{TopAngle}
    top_cg_pro_dihedrals::Vector{TopDihedral}
    top_cg_pro_aicg13::Vector{TopAngle}
    top_cg_pro_aicg14::Vector{TopDihedral}
    top_cg_pro_go_contact::Vector{TopContact}
    param_cg_pro_e_13::Vector{Float64}
    param_cg_pro_e_14::Vector{Float64}
    param_cg_pro_e_contact::Vector{Float64}
    # DNA
    top_cg_DNA_bonds::Vector{TopBond}
    top_cg_DNA_angles::Vector{TopAngle}
    top_cg_DNA_dih_Gaussian::Vector{TopDihedral}
    top_cg_DNA_dih_periodic::Vector{TopDihedral}
    param_cg_DNA_k_angles::Vector{Float64}
    # RNA
    top_cg_RNA_bonds::Vector{TopBond}
    top_cg_RNA_angles::Vector{TopAngle}
    top_cg_RNA_dihedrals::Vector{TopDihedral}
    top_cg_RNA_base_stack::Vector{TopContact}
    top_cg_RNA_base_pair::Vector{TopContact}
    top_cg_RNA_other_contact::Vector{TopContact}
    param_cg_RNA_k_bonds::Vector{Float64}
    param_cg_RNA_k_angles::Vector{Float64}
    param_cg_RNA_k_dihedrals::Vector{Float64}
    param_cg_RNA_e_base_stack::Vector{Float64}
    param_cg_RNA_e_base_pair::Vector{Float64}
    param_cg_RNA_e_other_contact::Vector{Float64}
    # protein-DNA
    top_cg_pro_DNA_pwmcos::Vector{TopPWMcos}
    # protein-RNA
    top_cg_pro_RNA_contact::Vector{TopContact}
    param_cg_pro_RNA_e_contact::Vector{Float64}
end

