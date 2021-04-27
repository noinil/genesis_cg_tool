###############################################################################
#                                                       _                     #
#                      __ _  _ __  ___    ___  _ __  __| |                    #
#                     / _` || '__|/ _ \  / __|| '__|/ _` |                    #
#                    | (_| || |  | (_) || (__ | |  | (_| |                    #
#                     \__, ||_|   \___/  \___||_|   \__,_|                    #
#                     |___/                                                   #
#                                                                             #
###############################################################################

using Printf

function write_grocrd(top::GenTopology, conf::Conformation, sys_name::AbstractString="", args::Dict{String, <:Any}=Dict{String, Any}())

    verbose  = get(args, "verbose", false)

    if length(sys_name) > 0
        system_name = sys_name
    else
        system_name = top.system_name
    end

    gro_name = system_name * ".gro"
    gro_file = open(gro_name, "w")

    cg_num_particles = conf.num_particle

    # measure box size
    (minx, miny, minz) = minimum(conf.coors, dims=2)
    (maxx, maxy, maxz) = maximum(conf.coors, dims=2)
    box_size = [maxx, maxy, maxz] - [minx, miny, minz]

    @printf(gro_file, "CG model %s, t = %16.3f \n", system_name, 0)
    @printf(gro_file, "%12d \n", cg_num_particles)

    for i_bead in 1 : cg_num_particles
        @printf(gro_file, "%5d%5s%5s%5d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n",
                top.top_atoms[i_bead].residue_index % 100000,
                top.top_atoms[i_bead].residue_name,
                top.top_atoms[i_bead].atom_name,
                i_bead % 100000,
                conf.coors[1 , i_bead] * 0.1,
                conf.coors[2 , i_bead] * 0.1,
                conf.coors[3 , i_bead] * 0.1,
                0.0, 0.0, 0.0)
    end
    @printf(gro_file, "%15.4f%15.4f%15.4f \n\n", box_size[1] * 0.1, box_size[2] * 0.1, box_size[3] * 0.1)

    close(gro_file)

    if verbose
        println(">           ... .gro: DONE!")
    end

end


function read_grocrd(grocrd_name::AbstractString)

    grocrd_lines = []
    for line in eachline(grocrd_name)
        push!(grocrd_lines, line)
    end

    num_particle = parse(Int, grocrd_lines[2])
    coors = zeros(Float64, (3, num_particle))
    for i in 1:num_particle
        words = split(grocrd_lines[i + 2][21:end])
        x = parse(Float64, words[1]) * 10.0
        y = parse(Float64, words[2]) * 10.0
        z = parse(Float64, words[3]) * 10.0
        coors[1, i] = x
        coors[2, i] = y
        coors[3, i] = z
    end
    conf = Conformation(num_particle, coors)

    return conf

end

