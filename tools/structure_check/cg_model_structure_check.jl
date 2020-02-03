#!/usr/bin/env julia

using Printf
using ArgParse

include("../../../src/lib/biomath.jl")
include("../../../src/lib/topology.jl")
include("../../../src/lib/constants.jl")
include("../../../src/lib/interactions.jl")
include("../../../src/lib/conformation.jl")
include("../../../src/lib/coarse_graining.jl")
include("../../../src/lib/parser_top.jl")

function main(top_filename::AbstractString)
    # ================================
    # Read in topology and coordinates
    # ================================

    mytop = read_grotop(top_filename)
    
end

main()

