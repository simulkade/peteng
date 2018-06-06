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

    # fluids, injection fluids, formation fluids
    fl = [:fluid_form, :fluid_inj]
    kw = ["formation_water", "injection_water"]
    for i in 1:2
        fluid = json_file["fluids"]
        if fluid["type"]=="constant"
            # extract the density, viscosity, and salinity data
            @eval $(fl[i])


        end

end