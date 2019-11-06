###############################################################################
#                             Molecule Structures                             #
###############################################################################

# ========
# All-atom
# ========

struct AAResidue
    name::String
    atoms::Array{Int64, 1}
end

struct AAChain
    id::Char
    segname::String
    moltype::Int
    residues::Array{Int64, 1}
end

struct AAMolecule
    atom_names::Array{String}
    atom_coors::Array{Float64, 2}
    residues::Array{AAResidue}
    chains::Array{AAChain}
end

struct AACGResidue
    res_idx::Int
    res_name::String
    atm_name::String
    atoms::Array{Int64, 1}
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

struct CGResidue
    name::String
    particles::Array{Int64, 1}
end

struct CGChain
    id::Char
    segname::String
    moltype::Int
    residues::Array{Int64, 1}
end

struct CGMolecule
    particle_names::Array{String}
    particle_coors::Array{Float64, 2}
    residues::Array{CGResidue}
    chains::Array{CGChain}
end

