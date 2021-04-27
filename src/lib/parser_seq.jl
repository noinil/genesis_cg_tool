###############################################################################
#     ____  _       _        __                            _   _              #
#    | __ )(_) ___ (_)_ __  / _| ___  _ __ _ __ ___   __ _| |_(_) ___ ___     #
#    |  _ \| |/ _ \| | '_ \| |_ / _ \| '__| '_ ` _ \ / _` | __| |/ __/ __|    #
#    | |_) | | (_) | | | | |  _| (_) | |  | | | | | | (_| | |_| | (__\__ \    #
#    |____/|_|\___/|_|_| |_|_|  \___/|_|  |_| |_| |_|\__,_|\__|_|\___|___/    #
#                                                                             #
###############################################################################

using Printf

function read_fasta(fasta_filename::String)
    mol_seqence = ""
    num_chain = 0
    seq_list = []
    for line in eachline(fasta_filename)
        if length(line) == 0
            continue
        end
        if line[1] == '>'
            if num_chain > 0 && length(mol_seqence) > 0
                push!(seq_list, mol_seqence)
            end
            mol_seqence = ""
            num_chain += 1
            continue
        end
        seq = strip(line)
        if length(seq) == 0
            continue
        end
        mol_seqence *= join(split(seq))
    end
    if num_chain > 0 && length(mol_seqence) > 0
        push!(seq_list, mol_seqence)
    end
    return (num_chain, seq_list[:])
end

function read_modified_pfm(pfm_filename::String)
    pfm = Dict()
    for line in eachline(pfm_filename)
        words = split(line)
        if length(words) < 1
            continue
        end
        w1 = words[1]
        if occursin(w1, "ACGT")
            local_list = []
            for p in words[2:end]
                push!( local_list, parse(Float64, p) )
            end
            pfm[w1] = local_list
        elseif in(w1, ["CHAIN_A", "CHAIN_B"])
            local_list = []
            for dna_id in words[2:end]
                push!( local_list, parse(Int, dna_id) )
            end
            pfm[w1] = local_list
        end
    end

    pfmat = [pfm["A"]  pfm["C"]  pfm["G"]  pfm["T"]]
    ppmat = pfmat ./ sum(pfmat, dims=2)
    pwmat0 = -log.(ppmat)
    pwmat = pwmat0 .- sum(pwmat0, dims=2) ./ 4

    return (pwmat, pfm["CHAIN_A"], pfm["CHAIN_B"])
end

function write_sequence(aa_molecule::AAMolecule, system_name::String)

    chain_id_set = "_ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"
    cg_seq_name = system_name * ".fasta"
    cg_seq_file = open(cg_seq_name, "w")

    aa_residues   = aa_molecule.residues
    aa_chains     = aa_molecule.chains

    aa_num_residue = length(aa_residues)
    aa_num_chain   = length(aa_chains)

    for i_chain in 1:aa_num_chain
        chain = aa_chains[i_chain]
        mol_type = chain.moltype
        @printf(cg_seq_file,
                "> Chain %s : %s \n",
                chain_id_set[mod(i_chain, 63) + 1],
                MOL_TYPE_LIST[mol_type])

        for i_res in chain.residues
            res_name = aa_residues[i_res].name
            print(cg_seq_file, RES_SHORTNAME_DICT[res_name])
        end

        print(cg_seq_file, "\n")
    end

    close(cg_seq_file)
    println(">           ... sequence output : DONE!")
 
end


