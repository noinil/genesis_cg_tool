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


# ==============
# Coarse-grained
# ==============

struct CGChain
    first::Int
    last::Int
    moltype::Int
    segname::String
end



