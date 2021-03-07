###############################################################################
#                             Molecule Structures                             #
###############################################################################

# ========
# All-atom
# ========

struct AAResidue
    name::String
    atoms::Vector{Int64}
end

struct AAChain
    id::String
    segname::String
    moltype::Int
    residues::Vector{Int64}
end

struct AAMolecule
    atom_names::Vector{String}
    atom_coors::Array{Float64, 2}
    residues::Vector{AAResidue}
    chains::Vector{AAChain}
end

struct AACGResidue
    res_idx::Int
    res_name::String
    atm_name::String
    atoms::Vector{Int64}
end

struct AACGChain
    first::Int
    last::Int
    moltype::Int
    segname::String
end

# ==============
# Coarse-grained
# ==============

# struct CGResidue
#     name::String
#     particles::Vector{Int64}
# end

# struct CGChain
#     id::Char
#     segname::String
#     moltype::Int
#     residues::Vector{Int64}
# end

# struct CGMolecule
#     particle_names::Vector{String}
#     particle_coors::Array{Float64, 2}
#     residues::Vector{CGResidue}
#     chains::Vector{CGChain}
# end

