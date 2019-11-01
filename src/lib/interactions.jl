module GCGInteractions

# --------------------
# AICG2+ Protein Model
# --------------------

function is_protein_backbone(atom_name::String)
    if in(atom_name, ("N", "C", "O", "OXT", "CA"))
        return true
    end
    return false
end

function is_protein_hb_donor(atom_name::String, res_name::String)
    if atom_name[1] == 'N'
        return true
    elseif atom_name[1] == 'S' && res_name == "CYS"
        return true
    elseif atom_name[1] == 'O'
        if  ( res_name == "SER" && atom_name == "OG"  ) ||
            ( res_name == "THR" && atom_name == "OG1" ) ||
            ( res_name == "TYR" && atom_name == "OH"  )
            return true
        end
    end
    return false
end

function is_protein_hb_acceptor(atom_name::String)
    if atom_name[1] == 'O' || atom_name[1] == 'S'
        return true
    end
    return false
end

function is_protein_cation(atom_name::String, res_name::String)
    if atom_name[1] == 'N'
        if  ( res_name == "ARG" && atom_name == "NH1" ) ||
            ( res_name == "ARG" && atom_name == "NH2" ) ||
            ( res_name == "LYS" && atom_name == "NZ"  )
            return true
        end
    end
    return false
end

function is_protein_anion(atom_name::String, res_name::String)
    if atom_name[1] == 'O'
        if  ( res_name == "GLU" && atom_name == "OE1" ) ||
            ( res_name == "GLU" && atom_name == "OE2" ) ||
            ( res_name == "ASP" && atom_name == "OD1" ) ||
            ( res_name == "ASP" && atom_name == "OD2" )
            return true
        end
    end
    return false
end

function is_protein_hb_pair(atom_name_1::String, res_name_1::String, atom_name_2::String, res_name_2::String)
    if  is_protein_hb_acceptor(atom_name_1) &&
        is_protein_hb_donor(atom_name_2, res_name_2)
        return true
    elseif is_protein_hb_acceptor(atom_name_2) &&
        is_protein_hb_donor(atom_name_1, res_name_1)
        return true
    end
    return false
end

function is_protein_sb_pair(atom_name_1::String, res_name_1::String, atom_name_2::String, res_name_2::String)
    if  is_protein_cation(atom_name_1, res_name_1) &&
        is_protein_anion(atom_name_2,  res_name_2)
        return true
    elseif is_protein_cation(atom_name_2, res_name_2) &&
        is_protein_anion(atom_name_1,  res_name_1)
        return true
    end
    return false
end

function is_protein_nonsb_charge_pair(atom_name_1::String, res_name_1::String, atom_name_2::String, res_name_2::String)
    if  is_protein_cation(atom_name_1, res_name_1) ||
        is_protein_anion(atom_name_1,  res_name_1) ||
        is_protein_cation(atom_name_2, res_name_2) ||
        is_protein_anion(atom_name_2,  res_name_2)
        return true
    end
    return false
end


end
