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
    cg_resid_name::Array{String}
    cg_resid_index::Array{Int}
    cg_bead_name::Array{String}
    cg_bead_type::Array{String}
    cg_bead_charge::Array{<:Real}
    cg_bead_mass::Array{<:Real}
    cg_chain_id::Array{Int}
    cg_seg_name::Array{String}
    # protein
    top_cg_pro_bonds::Array{TopBond}
    top_cg_pro_angles::Array{TopAngle}
    top_cg_pro_dihedrals::Array{TopDihedral}
    top_cg_pro_aicg13::Array{TopAngle}
    top_cg_pro_aicg14::Array{TopDihedral}
    top_cg_pro_go_contact::Array{TopContact}
    param_cg_pro_e_13::Array{Float64}
    param_cg_pro_e_14::Array{Float64}
    param_cg_pro_e_contact::Array{Float64}
    # DNA
    top_cg_DNA_bonds::Array{TopBond}
    top_cg_DNA_angles::Array{TopAngle}
    top_cg_DNA_dih_Gaussian::Array{TopDihedral}
    top_cg_DNA_dih_periodic::Array{TopDihedral}
    param_cg_DNA_k_angles::Array{Float64}
    # RNA
    top_cg_RNA_bonds::Array{TopBond}
    top_cg_RNA_angles::Array{TopAngle}
    top_cg_RNA_dihedrals::Array{TopDihedral}
    top_cg_RNA_base_stack::Array{TopContact}
    top_cg_RNA_base_pair::Array{TopContact}
    top_cg_RNA_other_contact::Array{TopContact}
    param_cg_RNA_k_bonds::Array{Float64}
    param_cg_RNA_k_angles::Array{Float64}
    param_cg_RNA_k_dihedrals::Array{Float64}
    param_cg_RNA_e_base_stack::Array{Float64}
    param_cg_RNA_e_base_pair::Array{Float64}
    param_cg_RNA_e_other_contact::Array{Float64}
    # protein-DNA
    top_cg_pro_DNA_pwmcos::Array{TopPWMcos}
    # protein-RNA
    top_cg_pro_RNA_contact::Array{TopContact}
    param_cg_pro_RNA_e_contact::Array{Float64}
end

