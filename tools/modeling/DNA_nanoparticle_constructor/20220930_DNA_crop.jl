#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/gcj.jl")
using ArgParse

function DNA_crop(DNA_name, nt1, nt2, out_name)
    # -------------
    # index mapping
    # -------------
    n_bp_old = 21

    n_DNA_strand1_nt = nt1 > n_bp_old ? n_bp_old : nt1
    n_DNA_strand2_nt = nt2 > n_bp_old ? n_bp_old : nt2

    # -----------------
    # top and crd files
    # -----------------
    tmp_top_name = @sprintf("%s.top", DNA_name)
    tmp_crd_name = @sprintf("%s.gro", DNA_name)
    mytop = read_grotop(tmp_top_name)
    mycrd = read_grocrd(tmp_crd_name)




    # =============================
    # Make a new topology structure
    # =============================

    n_DNA_strand1_atom = n_DNA_strand1_nt > 0 ? n_DNA_strand1_nt * 3 - 1 : 0
    n_DNA_strand2_atom = n_DNA_strand2_nt > 0 ? n_DNA_strand2_nt * 3 - 1 : 0
    n_atom_DNA_new = n_DNA_strand1_atom + n_DNA_strand2_atom

    old_indices = vcat([1:n_DNA_strand1_atom...], [n_bp_old * 3:n_bp_old * 3 - 1 + n_DNA_strand2_atom...])
    # new_indices = [1:n_atom_DNA_new...]

    top_default_params             = GenTopDefault(0, 0, false, 0.0, 0.0)
    top_default_atomtype           = Vector{GenTopAtomType}(undef, 0)
    top_default_CGDNA_bp           = Vector{GenTopCGDNABasepairType}(undef, 0)
    top_default_CGDNA_bs           = Vector{GenTopCGDNABasestackType}(undef, 0)
    top_default_CGDNA_cs           = Vector{GenTopCGDNABasecrossType}(undef, 0)
    top_default_CGDNA_exv          = Vector{GenTopCGDNAExvType}(undef, 0)
    top_default_CGPro_flx_angle    = Vector{GenTopCGProAICGFlexAngleType}(undef, 0)
    top_default_CGPro_flx_dihedral = Vector{GenTopCGProAICGFlexDihedralType}(undef, 0)

    global_index_2_local_index     = Vector{Int}(undef, 0)
    global_index_2_local_molid     = Vector{Int}(undef, 0)
    top_atoms                      = Vector{GenTopAtom}(undef, 0)
    top_bonds                      = Vector{GenTopBond}(undef, 0)
    top_angles                     = Vector{GenTopAngle}(undef, 0)
    top_dihedrals                  = Vector{GenTopDihedral}(undef, 0)
    top_pairs                      = Vector{GenTopPair}(undef, 0)
    top_exclusions                 = Vector{GenTopExclusion}(undef, 0)
    top_pwmcos                     = Vector{GenTopPWMcos}(undef, 0)
    top_pwmcosns                   = Vector{GenTopPWMcos}(undef, 0)
    top_idr_hps                    = Vector{GenTopRegion}(undef, 0)
    top_idr_kh                     = Vector{GenTopRegion}(undef, 0)
    top_mol_list                   = Vector{GenTopMolList}(undef, 0)

    println(" > Preparing topology...")
    # ---------
    # [ atoms ]
    # ---------
    for (k, k_bead) in enumerate(old_indices)
        # new_atom = GenTopAtom(i_bead, a_type, r_indx, r_name, a_name, f_type, charge, mass, c_id, s_name)
        my_atom = mytop.top_atoms[k_bead]
        n_dna_res = my_atom.residue_index
        new_atom = GenTopAtom(k, my_atom.atom_type,
                              n_dna_res,
                              my_atom.residue_name,
                              my_atom.atom_name,
                              AICG_ATOM_FUNC_NR,
                              my_atom.charge,
                              my_atom.mass,
                              my_atom.chain_id, my_atom.seg_name)
        push!(top_atoms, new_atom)
    end
    # ---------
    # [ bonds ]
    # ---------
    n_bond_old = length(mytop.top_bonds)
    for i_bond in 1:n_bond_old
        b = mytop.top_bonds[i_bond]
        if !( b.i in old_indices ) || !(b.j in old_indices)
            continue
        end
        new_bond = GenTopBond(indexin(b.i, old_indices)[1], indexin(b.j, old_indices)[1], b.function_type, b.r0, b.coef)
        push!(top_bonds, new_bond)
    end

    # ----------
    # [ angles ]
    # ----------
    n_angle_old = length(mytop.top_angles)
    for i_angle in 1:n_angle_old
        a = mytop.top_angles[i_angle]
        if !( a.i in old_indices ) || !(a.j in old_indices) || !(a.k in old_indices)
            continue
        end
        new_angle = GenTopAngle(indexin(a.i, old_indices)[1], indexin(a.j, old_indices)[1], indexin(a.k, old_indices)[1], DNA3SPN_ANG_FUNC_TYPE, a.a0, a.coef, 0.0)
        push!(top_angles, new_angle)
    end

    # -------------
    # [ dihedrals ]
    # -------------
    n_dihedral_old = length(mytop.top_dihedrals)
    for i_dihedral in 1:n_dihedral_old
        d = mytop.top_dihedrals[i_dihedral]
        if !( d.i in old_indices ) || !(d.j in old_indices) || !(d.k in old_indices) || !(d.l in old_indices)
            continue
        end
        new_dihedral = GenTopDihedral(indexin(d.i, old_indices)[1],
                                      indexin(d.j, old_indices)[1],
                                      indexin(d.k, old_indices)[1],
                                      indexin(d.l, old_indices)[1],
                                      d.function_type,
                                      d.d0, d.coef, d.w, d.n)
        push!(top_dihedrals, new_dihedral)
    end

    # -----------
    # Coordinates
    # -----------
    atom_coors = zeros(Float64, (3, n_atom_DNA_new))
    atom_coors[:, :] = mycrd.coors[:, old_indices]



    # ==========
    # new output
    # ==========
    println(" ------------------------------------------------------------ ")
    println(" > Output files...")
    if length(out_name) > 0
        mol_name = out_name
    else
        mol_name = @sprintf("new_%s", DNA_name)
    end
    mytop = GenTopology(mol_name, n_atom_DNA_new,
                        top_default_params,
                        top_default_atomtype,
                        top_default_CGDNA_bp,
                        top_default_CGDNA_bs,
                        top_default_CGDNA_cs,
                        top_default_CGDNA_exv,
                        top_default_CGPro_flx_angle,
                        top_default_CGPro_flx_dihedral,
                        global_index_2_local_index,
                        global_index_2_local_molid,
                        top_atoms,
                        top_bonds,
                        top_angles,
                        top_dihedrals,
                        top_pairs,
                        top_exclusions,
                        top_pwmcos,
                        top_pwmcosns,
                        top_idr_hps,
                        top_idr_kh,
                        top_mol_list)
    myconf = Conformation(n_atom_DNA_new, atom_coors)

    # -----------------
    # output some files
    # -----------------
    write_grotop(mytop, mol_name)
    write_grocrd(mytop, myconf, mol_name)
    write_mmCIF(mytop, myconf, mol_name)

end

if abspath(PROGRAM_FILE) == @__FILE__

    my_args = ArgParseSettings()

    @add_arg_table my_args begin
        "--DNA_file_name", "-D"
        help     = "File name for DNA topology and coordinates (basename of .top, .gro)"
        arg_type = String
        default  = ""

        "--output_file_name", "-o"
        help     = "File name for output"
        arg_type = String
        default  = ""

        "--nt1"
        help     = "Number of nucleotides in 1st strand"
        arg_type = Int
        default  = 1000000

        "--nt2"
        help     = "Number of nucleotides in 2nd strand"
        arg_type = Int
        default  = 1000000

        "--debug"
        help = "Debug mode"
        action = :store_true
    end

    main_args = parse_args(my_args)

    DNA_name = main_args["DNA_file_name"]
    out_name = main_args["output_file_name"]
    nt1 = main_args["nt1"]
    nt2 = main_args["nt2"]

    if nt1 < 0 || nt2 < 0
        println(" Please choose nt1 >= 0 and nt2 >= 0!")
        exit()
    end

    DNA_crop(DNA_name, nt1, nt2, out_name)
end
