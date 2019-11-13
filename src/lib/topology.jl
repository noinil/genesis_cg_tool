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

# -------------
# Default types
# -------------

struct GenTopDefault
    nonbonded_function_type::Int
    nonbonded_combination_rule::Int
    gen_pairs::Bool
    fudge_lj::Float64
    fudge_qq::Float64
end

struct GenTopAtomType
    name::String
    mass::Float64
    charge::Float64
    rmin::Float64
    eps::Float64
    n::Int
    ptype::String
end

struct GenTopCGDNABasestackType
    base_5::String
    base_3::String
    function_type::Int
    epsilon::Float64
    sigma::Float64
    theta::Float64
end

struct GenTopCGDNABasecrossType
    base_a::String
    base_b::String
    function_type::Int
    epsilon::Float64
    sigma::Float64
    theta::Float64
end

struct GenTopCGDNABasepairType
    base_a::String
    base_b::String
    function_type::Int
    theta1::Float64
    theta2::Float64
    theta3::Float64
    phi1::Float64
    sigma::Float64
    epsilon::Float64
end

struct GenTopCGDNAExvType
    base::String
    function_type::Int
    sigma::Float64
end

struct GenTopCGProAICGFlexAngleType
    # TODO
end

struct GenTopCGProAICGFlexDihedralType
    # TODO
end

# -----------------------
# Molecule specific types
# -----------------------

struct GenTopAtom
    atom_index::Int
    atom_type::String
    residue_index::Int
    residue_type::String
    atom_name::String
    function_type::Int
    charge::Float64
    mass::Float64
end

struct GenTopBond
    i::Int
    j::Int
    function_type::Int
    r0::Float64
    coef::Float64
end

struct GenTopAngle
    i::Int
    j::Int
    k::Int
    function_type::Int
    a0::Float64
    coef::Float64
    w::Float64
end

struct GenTopDihedral
    i::Int
    j::Int
    k::Int
    l::Int
    function_type::Int
    d0::Float64
    coef::Float64
    w::Float64
    n::Int
end

struct GenTopPair
    i::Int
    j::Int
    function_type::Int
    r0::Float64
    coef::Float64
end

struct GenTopExclusion
    i::Int
    j::Int
end

struct GenTopPWMcos
    i::Int
    function_type::Int
    r0::Float64
    theta1::Float64
    theta2::Float64
    theta3::Float64
    ene_A::Float64
    ene_C::Float64
    ene_G::Float64
    ene_T::Float64
    gamma::Float64
    eps::Float64
end

struct GenTopMolecule
    mol_name::String
    nonlocal_interval::Int
    num_atom::Int

    top_atoms::Vector{GenTopAtom}
    top_bonds::Vector{GenTopBond}
    top_angles::Vector{GenTopAngle}
    top_dihedrals::Vector{GenTopDihedral}
    top_pairs::Vector{GenTopPair}
    top_exclusions::Vector{GenTopExclusion}
    top_pwmcos::Vector{GenTopPWMcos}
end

struct GenTopMolList
    mol_name::String
    count::Int
end


struct GenTopology
    system_name::String
    num_atom::Int

    top_default_params::GenTopDefault
    top_default_atomtype::Vector{GenTopAtomType}
    top_default_CGDNA_bp::Vector{GenTopCGDNABasepairType}
    top_default_CGDNA_bs::Vector{GenTopCGDNABasestackType}
    top_default_CGDNA_cs::Vector{GenTopCGDNABasecrossType}
    top_default_CGDNA_exv::Vector{GenTopCGDNAExvType}
    top_default_CGPro_flx_angle::Vector{GenTopCGProAICGFlexAngleType}
    top_default_CGPro_flx_dihedral::Vector{GenTopCGProAICGFlexDihedralType}

    global_index_2_local_index::Vector{Int}
    global_index_2_local_molid::Vector{Int}
    top_atoms::Vector{GenTopAtom}
    top_bonds::Vector{GenTopBond}
    top_angles::Vector{GenTopAngle}
    top_dihedrals::Vector{GenTopDihedral}
    top_pairs::Vector{GenTopPair}
    top_exclusions::Vector{GenTopExclusion}
    top_pwmcos::Vector{GenTopPWMcos}

    top_mol_list::Vector{GenTopMolList}
   
end

