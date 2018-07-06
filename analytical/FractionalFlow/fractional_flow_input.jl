# Read the input json file and assign it to the FractionalFlow
# structures

"""
read_frac_flow(input_file::AbstractString)
This function reads the fractional flow input in JSON
format and returns a set of structures that defines 
an injection problem
"""
function read_frac_flow_input(input_file::AbstractString)
    json_file = JSON.parsefile(input_file)

    # Domain (i.e., the core)
    L_core  = json_file["core"]["length"]
    D_core  = json_file["core"]["diameter"]
    k       = json_file["core"]["permeability"]
    φ       = json_file["core"]["porosity"]
    a_core  = json_file["core"]["surface_area"]
    ρ_core  = json_file["core"]["density"]
    core_props = CoreProperties(L_core, D_core, φ, k, a_core, ρ_core)

    # relative permeabilities (injected, formation)
    # slightly fancy with @eval
    rp = [:rel_perm_form, :rel_perm_inj]
    kw = ["formation", "injection"]
    for i in 1:2
        rel_perm = json_file["core"]["relative_permeability"][kw[i]]
        if rel_perm["type"]=="Corey"
            kro0 = rel_perm["kro0"]
            krw0 = rel_perm["krw0"]
            no = rel_perm["no"]
            nw = rel_perm["nw"]
            swc = rel_perm["swc"]
            sor = rel_perm["sor"]
            @eval $(rp[i]) = CoreyRelativePermeability(krw0, kro0, swc, sor, nw, no)
        # elseif [other rel perm types for the future]
        end
    end

    # fluids, injection fluids, formation fluids
    
    

end

function read_brine(brine::Dict{String, Any}, brine_name::AbstractString="brine", 
    brine_number::Integer=1)
    
    brine_sal = Array(Brine, length(brine["salinity"]))
    i = 0
    for k in keys(brine["salinity"])
        i+=1
        brine_sal[i].species = k
        brine_sal[i].concentration = brine["salinity"][k]
    end
    return Brine(brine_name, brine_number, 
        get(brine, "temperature", 300.0),
        get(brine, "pressure", 1.01e5),
        get(brine, "unit", "mol/kgw"),
        brine_sal,
        get(brine, "pH", 7.0),
        get(brine, "density", 1000.0),
        get(brine, "viscosity", 0.001)
    )
end
