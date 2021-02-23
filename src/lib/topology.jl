###############################################################################
#                                   Topology                                  #
###############################################################################

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
    residue_name::String
    atom_name::String
    function_type::Int
    charge::Float64
    mass::Float64
    chain_id::Int
    seg_name::String
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

struct GenTopRegion
    istart::Int
    iend::Int
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
    top_pwmcosns::Vector{GenTopPWMcos}
    top_idr_hps::Vector{GenTopRegion}
    top_idr_kh::Vector{GenTopRegion}
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
    top_pwmcosns::Vector{GenTopPWMcos}
    top_idr_hps::Vector{GenTopRegion}
    top_idr_kh::Vector{GenTopRegion}

    top_mol_list::Vector{GenTopMolList}
   
end

