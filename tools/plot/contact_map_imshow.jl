#!/usr/bin/env julia

using ArgParse
using PyPlot

function plot_contact_map(itp_filename, do_dist_matrix, cmap_name, do_debug)
    i_step = 1
    println("============================================================")
    println("> Step $(i_step): read [ atoms ].")
    # -----------------------
    # read in number of atoms
    # -----------------------
    num_atom = 0
    do_atom_read = false
    for line in eachline(itp_filename)
        if do_atom_read
            if length(line) == 0
                break
            end
            if line[1] == ';'
                continue
            else
                words = split(line)
                if length(words) == 8
                    num_atom += 1
                elseif line[1] == '['
                    break
                end
            end
        end
        if startswith(line, "[ atoms ]")
            do_atom_read = true
        end
    end

    println("------------------------------------------------------------")
    println("          In total: $(num_atom) CG particles.")

    i_step += 1
    println("============================================================")
    println("> Step $(i_step): read [ pairs ].")
    # -------------------
    # read in contact map
    # -------------------
    num_pair = 0
    contact_matrix = zeros(Float64, ( num_atom, num_atom ))
    r_min = 1e50
    r_max = -1e50
    do_pair_read = false
    for line in eachline(itp_filename)
        if do_pair_read
            if length(line) == 0
                break
            end
            if line[1] == ';'
                continue
            else
                words = split(line)
                if length(words) == 5
                    num_pair += 1
                    i = parse(Int, words[1])
                    j = parse(Int, words[2])
                    r = parse(Float64, words[4])
                    contact_matrix[i, j] = do_dist_matrix ? r : 1.0
                    contact_matrix[j, i] = contact_matrix[i, j]
                    r_min = r < r_min ? r : r_min
                    r_max = r > r_max ? r : r_max
                elseif line[1] == '['
                    break
                end
            end
        end
        if startswith(line, "[ pairs ]")
            do_pair_read = true
        end
    end
    println("------------------------------------------------------------")
    println("          In total: $(num_pair) contacts.")

    i_step += 1
    println("============================================================")
    println("> Step $(i_step): plotting...")
    # -----------
    # Plotting...
    # -----------
    fig, ax = plt.subplots(figsize=(16, 16))
    if do_dist_matrix
        ax.imshow(contact_matrix, cmap=cmap_name, vmin=r_min, vmax=r_max, aspect="equal", origin="lower")
    else
        ax.imshow(contact_matrix, cmap=cmap_name, aspect="equal", origin="lower")
    end

    svg_filename = itp_filename[1:end-4] * "_contact_map.svg"
    savefig(svg_filename)
    println("------------------------------------------------------------")
    println("          Saved to: $(svg_filename).")
    println("[1;32m DONE! [0m ")
    println("============================================================")

end


# =============================
# Parsing Commandline Arguments
# =============================
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "itp"
        help     = "Gromacs style .itp file name."
        required = true
        arg_type = String

        "--distance-matrix"
        help     = "Plot with distance based color-map."
        action   = :store_true

        "--cmap"
        help     = "Color-map."
        default  = "gray"
        arg_type = String

        "--debug"
        help     = "DEBUG."
        action   = :store_true
    end

    return parse_args(s)
end

# ====
# Main
# ====
function main()
    args = parse_commandline()

    plot_contact_map(args["itp"],
                     args["distance-matrix"],
                     args["cmap"],
                     args["debug"])

end

main()
