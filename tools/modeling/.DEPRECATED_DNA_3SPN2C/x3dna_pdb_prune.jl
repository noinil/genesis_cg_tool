#!/usr/bin/env julia

using ArgParse

# =============================
# Parsing Commandline Arguments
# =============================
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin

        "pdb"
        help     = "x3dna DNA pdb file name."
        required = true
        arg_type = String

    end

    return parse_args(s)
end

# ====
# Main
# ====
function main()

    args = parse_commandline()

    old_pdb_name = args["pdb"]

    new_pdb_name = old_pdb_name[1:end-4] * "_new.pdb"

    out_file = open(new_pdb_name, "w")

    tmp_id_chain = ' '
    for line in eachline(args["pdb"])
        if startswith(line, "ATOM")
            chain_id          = line[22]
            if chain_id != tmp_id_chain
                if tmp_id_chain != ' '
                    println(out_file, "TER")
                end
                tmp_id_chain = chain_id
            end
            println(out_file, line)
        elseif startswith(line, "TER") || startswith(line, "END")
            println(out_file, line)
        end
    end

    close(out_file)

end

main()
