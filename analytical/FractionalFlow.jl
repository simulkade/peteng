module FractionalFlow

using Roots, Dierckx, PyPlot, JFVM
include("CoreyFunctions.jl")
include("water_flood_fvm.jl")
include("fractional_flow_cases.jl")
CF = CoreyFunctions

# types:
struct CoreyRelativePermeability
    krw0::Real
    kro0::Real
    swc::Real
    sor::Real
    nw::Real
    no::Real
end

struct CoreProperties
    length::Real
    diameter::Real
    porosity::Real
    permeability::Real
end

struct CoreFlooding
    injection_velocity::Real # m/s
    injected_pore_volume::Real
    back_pressure::Real
    initial_water_saturation::Real
    injected_water_saturation::Real
end

struct Fluids
    oil_viscosity::Real
    water_viscosity::Real
end

"""
A structure that stores the fractional flow curves and the final solution
"""
struct FracFlowResults
    fractional_flow_functions::Array{Function, 1}
    shock_lines::Array{Real, 2}
    recovery_pv::Array{Real, 2}
    recovery_time::Array{Real, 2}
    saturation_profile_xt::Array{Real, 2}
end

# functions

"""
"""
function oil_water_rel_perms(;krw0=0.4, kro0=0.9, 
    swc=0.15, sor=0.2, nw=2.0, no = 2.0)
    return CoreyRelativePermeability(krw0, kro0, swc, sor, nw, no)
end

function oil_water_fluids(;mu_water=1e-3, mu_oil=2e-3)
    return Fluids(mu_oil, mu_water)
end

function core_flooding(;u_inj=1.15e-5, pv_inject=5.0, p_back=1e5, sw_init=0.2, sw_inj=1.0, rel_perms=[])
    swc = sw_init
    if typeof(rel_perms)==CoreyRelativePermeability
        if rel_perms.swc>sw_init
            swc = rel_perms.swc
        end
    end
    return CoreFlooding(u_inj, pv_inject, p_back, swc, sw_inj)
end

function core_properties(;L=0.15, D=0.03, φ=0.3, k=1e-12)
    return CoreProperties(L, D, φ, k)
end



function rel_perm_functions(rel_perm::CoreyRelativePermeability)
    kro0 = rel_perm.kro0
    krw0 = rel_perm.krw0
    sor  = rel_perm.sor
    swc = rel_perm.swc
    no = rel_perm.no
    nw = rel_perm.nw

    kro = sw -> CF.kro(sw, kro0, sor, swc, no)
    krw = sw -> CF.krw(sw, krw0, sor, swc, nw)
    dkrwdsw = sw -> CF.dkrwdsw(sw, krw0, sor, swc, nw)
    dkrodsw = sw -> CF.dkrodsw(sw, kro0, sor, swc, no)
  
    return krw, kro, dkrwdsw, dkrodsw
end

function fractional_flow_function(rel_perms, fluids)
    krw, kro, dkrwdsw, dkrodsw = rel_perm_functions(rel_perms)
    muw = fluids.water_viscosity
    muo = fluids.oil_viscosity

    fw = sw -> ((krw(sw)/muw)/(krw(sw)/muw+kro(sw)/muo))
    dfw = sw -> ((dkrwdsw(sw)/muw*(krw(sw)/muw+kro(sw)/muo)-
    (dkrwdsw(sw)/muw+dkrodsw(sw)/muo)*krw(sw)/muw)/
    (kro(sw)/muo+krw(sw)/muw)^2)
    
    return fw, dfw
end

function tangent_line_saturation(rel_perms, fluids, sw_fw_point)
    swc = rel_perms.swc
    sor = rel_perms.sor
    sw0 = sw_fw_point[1]
    fw0 = sw_fw_point[2]
    fw, dfw = fractional_flow_function(rel_perms, fluids)

    f_shock = sw -> (dfw(sw)-(fw(sw)-fw0)/(sw-sw0))
    sw_tmp = collect(linspace(swc+eps(), 1-sor-eps(), 1000))
    # ind_min = indmin(abs.(f_shock.(sw_tmp)))
    ind_max = indmax(dfw.(sw_tmp))
    sw_max = sw_tmp[ind_max]
    ind_min = indmin(abs.(f_shock.(sw_tmp[ind_max:end])))
    sw_shock = sw_tmp[ind_min]
    try
        if f_shock(sw_max)*f_shock(1-sor-eps())>0
            sw_shock = fzero(f_shock, sw_shock) #  ftol = 1e-10, xtol = 1e-7)
        else
            sw_shock = fzero(f_shock, [sw_max,1-sor-eps()])
        end
    catch()
        info("shock front saturation is estimated: $sw_shock, error is $(f_shock(sw_shock))")
    end
end

function cross_point_saturation(fw, rel_perms, point1, point2)
    # create the equation for line
    sor = rel_perms.sor
    swc = rel_perms.swc
    sw_tmp = collect(linspace(swc+eps(), 1-sor-eps(), 1000))
    m = (point2[2]-point1[2])/(point2[1]-point1[1])
    b = point1[2]-m*point1[1]
    f = sw -> m*sw+b-fw(sw)
    f_tmp = f.(sw_tmp)
    sw_neg = sw_tmp[indmin(f_tmp)]
    sw_pos = sw_tmp[indmax(f_tmp)]
    sw_left, sw_right = minmax(sw_neg, sw_pos)
    sw_cross = sw_tmp[indmin(abs.(f_tmp))]
    try
        if f(sw_left)*f(sw_right)>0
            sw_cross = fzero(f, sw_cross)
        else
            sw_cross = fzero(f, [sw_left, sw_right])
        end
    catch()
        info("Cross point is estimated: $sw_cross, error is $(f(sw_cross))")
    end
    return sw_cross
end

function water_flood(core_props, fluids, rel_perms, core_flood)
    fw, dfw = fractional_flow_function(rel_perms, fluids)
    sw_init = core_flood.initial_water_saturation
    sw_shock = tangent_line_saturation(rel_perms, fluids, (sw_init, fw(sw_init)))
    println("sw_shock = $sw_shock")
    t_D_BT = 1/dfw(sw_shock) # breakthrough (BT) time [#PV]
    println("breakthrough time = $t_D_BT")
    
    # construct the recovery factor curve versus the # of PV
    R = zeros(1)
    pv_R = zeros(1)
    # at breakthrough
    push!(R, (1-fw(sw_init))*t_D_BT/(1-sw_init)) # recovery at BT
    push!(pv_R, t_D_BT) # BT time
    # after breakthrough
    pv_inj = max(core_flood.injected_pore_volume, 2.0) # at least inject 2 pv
    f_sw = sw -> (pv_inj-1/dfw(sw)) # find the outlet saturation at the end of injection
    sw_max = fzero(f_sw, sw_shock)
    sw_tmp = linspace(sw_shock, sw_max, 100)
    t_D_tmp = 1./dfw.(sw_tmp)
    
    s_av_tmp = sw_tmp-(fw.(sw_tmp)-1).*t_D_tmp
    R_tmp = (s_av_tmp-sw_init)/(1-sw_init)
    append!(R, R_tmp)
    append!(pv_R, t_D_tmp)

    # saturation profile
    sor = rel_perms.sor
    sw_inj = core_flood.injected_water_saturation
    phi = core_props.porosity
    ut = core_flood.injection_velocity
    s1 = collect(linspace(min(sw_inj, 1-sor-eps()), sw_shock-eps(), 100))
    xt_s1 = dfw.(s1)
    xt_shock = dfw(sw_shock)
    xt_prf=[xt_s1; xt_shock; xt_shock+eps(); 2*xt_shock]
    sw_prf=[s1; sw_shock; sw_init; sw_init]

    return pv_R, R, xt_prf, sw_prf # for the time being to test the code

end

"""
Low salinity water injection into a high salinity aquifer.
It is assumed that the low salinity fluid has a different relative permeability
and viscosity. Therefore, two fluids objects needs to be defined.
In this model, it is assumed that the core is not initially fooded with the 
reservoir brine, i.e., low salinity is the secondary flooding.
"""
function low_sal_water_flood(core_props, fluids_ls, fluids_hs, rel_perms_hs, 
    rel_perms_ls, core_flood)
    pv_R, R, xt_prf, sw_prf= single_ion_adsorption_water_flood(core_props, fluids_ls, fluids_hs, 
    rel_perms_hs, rel_perms_ls, core_flood, 0.0)
    return pv_R, R, xt_prf, sw_prf # for the time being to test the code
end

function single_ion_adsorption_water_flood(core_props, fluids_ls, fluids_hs, rel_perms_hs, 
    rel_perms_ls, core_flood, eq_const)
    # construct the fractional flow curves
    fw_ls, dfw_ls = fractional_flow_function(rel_perms_ls, fluids_ls)
    fw_hs, dfw_hs = fractional_flow_function(rel_perms_hs, fluids_hs)
    # low sal shock (tangent line from (0,0))
    sw_shock_ls = tangent_line_saturation(rel_perms_ls, fluids_ls, (-eq_const, 0.0))
    println("low sal sw_shock = $sw_shock_ls")
    t_D_BT_ls = 1/dfw_ls(sw_shock_ls) # breakthrough (BT) time [#PV]
    println("low sal breakthrough time = $t_D_BT_ls")
    # High sal shock (cross point between the ls tangent and the hs fw)
    sw_shock_hs = cross_point_saturation(fw_hs, rel_perms_hs, (-eq_const, 0.0), (sw_shock_ls, fw_ls(sw_shock_ls)))
    println("high sal sw_shock = $sw_shock_hs")
    sw_init = core_flood.initial_water_saturation    
    t_D_BT_hs = (sw_shock_hs-sw_init)/(fw_hs(sw_shock_hs)-fw_hs(sw_init)) # breakthrough (BT) time [#PV]
    println("high sal breakthrough time = $t_D_BT_hs")
    
    # construct the recovery factor curve versus the # of PV
    R = zeros(1)
    pv_R = zeros(1)
    # at breakthrough of the hs brine
    push!(R, (1-fw_hs(sw_init))*t_D_BT_hs/(1-sw_init)) # recovery at BT
    push!(pv_R, t_D_BT_hs) # BT time
    # at breakthrough of the ls brine
    push!(R, R[end]+(1-fw_hs(sw_shock_hs))*(t_D_BT_ls-t_D_BT_hs)/(1-sw_init)) # recovery at BT
    push!(pv_R, t_D_BT_ls) # BT time

    # after breakthrough
    pv_inj = max(core_flood.injected_pore_volume, 2.0) # at least inject 2 pv
    f_sw = sw -> (pv_inj-1/dfw_ls(sw)) # find the outlet saturation at the end of injection
    sw_max = fzero(f_sw, sw_shock_ls)
    sw_tmp = linspace(sw_shock_ls, sw_max, 100)
    t_D_tmp = 1./dfw_ls.(sw_tmp)
    
    s_av_tmp = sw_tmp-(fw_ls.(sw_tmp)-1).*t_D_tmp
    R_tmp = (s_av_tmp-sw_init)/(1-sw_init)
    append!(R, R_tmp)
    append!(pv_R, t_D_tmp)

    # saturation profile
    sor = rel_perms_ls.sor
    sw_inj = core_flood.injected_water_saturation
    phi = core_props.porosity
    ut = core_flood.injection_velocity
    s1 = collect(linspace(min(sw_inj, 1-sor-eps()), sw_shock_ls-eps(), 100))
    xt_s1 = dfw_ls.(s1)
    xt_shock_ls = 1/t_D_BT_ls
    xt_shock_hs = 1/t_D_BT_hs
    xt_prf=[xt_s1; xt_shock_ls; xt_shock_ls+eps(); xt_shock_hs; xt_shock_hs+eps(); 1/0.3]
    sw_prf=[s1; sw_shock_ls; sw_shock_hs; sw_shock_hs; sw_init; sw_init]

    return pv_R, R, xt_prf, sw_prf # for the time being to test the code
end

function single_ion_adsorption_tertiary_water_flood(core_props, fluids_ls, fluids_hs, rel_perms_hs, 
    rel_perms_ls, core_flood, eq_const)
    # ====================================================================
    # This function is not done yet. The trick is to draw a tangent from 
    # the high sal 1-sor_hs to the low sal curve to find the shock front 
    # saturation. The low sal shock speed is calculated at this saturation
    # by fw/(sw+eq_const)
    # The implementation needs some if statements to make sure that the 
    # shock front saturation calculated by the tangent line from -eq_const
    # to the low sal fw is indeed lower than 1-sor_hs
    # ====================================================================
    # construct the fractional flow curves
    fw_ls, dfw_ls = fractional_flow_function(rel_perms_ls, fluids_ls)
    fw_hs, dfw_hs = fractional_flow_function(rel_perms_hs, fluids_hs)
    # low sal shock (tangent line from (0,0))
    sw_shock_ls = tangent_line_saturation(rel_perms_ls, fluids_ls, (-eq_const, 0.0))
    println("low sal sw_shock = $sw_shock_ls")
    t_D_BT_ls = 1/dfw_ls(sw_shock_ls) # breakthrough (BT) time [#PV]
    println("low sal breakthrough time = $t_D_BT_ls")
    # This is speculation: sometimes the calculated low sal shock sw is lower than 1-sor_ls
    # which is impossible. Therefore, I correct it here, and also correct the low sal
    # breakthrough time. Therefore, the following correction is applied:
    sor_hs = rel_perms_hs.sor
    if (sw_shock_ls<1-sor_hs)
        info("the calculated shock front $sw_shock_ls is less than $(1-sor_hs). Correection is needed")
        sw_shock_ls = tangent_line_saturation(rel_perms_ls, fluids_ls, (1-sor_hs, 0.0))
        println("Corrected: low sal sw_shock = $sw_shock_ls")
        t_D_BT_ls = dfw_ls(sw_shock_ls)
        info("low sal breakthrough time corrected to $t_D_BT_ls_w")
    end
    
    # High sal shock (cross point between the ls tangent and the hs fw)
    sw_shock_hs = cross_point_saturation(fw_hs, rel_perms_hs, (-eq_const, 0.0), (sw_shock_ls, fw_ls(sw_shock_ls)))
    println("high sal sw_shock = $sw_shock_hs")
    t_D_BT_hs = (sw_shock_hs-(1-sor_hs))/(fw_hs(sw_shock_hs)-fw_hs(1-sor_hs)) # breakthrough (BT) time [#PV]
    println("high sal breakthrough time = $t_D_BT_hs")
    
    

    # construct the recovery factor curve versus the # of PV
    R = zeros(1)
    pv_R = zeros(1)
    # at breakthrough of the hs brine
    sw_init = core_flood.initial_water_saturation
    push!(R, 0.0) # recovery at BT
    push!(pv_R, t_D_BT_hs) # BT time
    # at breakthrough of the ls brine
    push!(R, (1-fw_hs(sw_shock_hs))*(t_D_BT_ls-t_D_BT_hs)/(1-sw_init)) # recovery at BT
    push!(pv_R, t_D_BT_ls) # BT time

    # after breakthrough
    pv_inj = max(core_flood.injected_pore_volume, 2.0) # at least inject 2 pv
    f_sw = sw -> (pv_inj-1/dfw_ls(sw)) # find the outlet saturation at the end of injection
    sw_max = fzero(f_sw, sw_shock_ls)
    sw_tmp = linspace(sw_shock_ls, sw_max, 100)
    t_D_tmp = 1./dfw_ls.(sw_tmp)
    
    s_av_tmp = sw_tmp-(fw_ls.(sw_tmp)-1).*t_D_tmp
    R_tmp = (s_av_tmp-(1-sor_hs))/(1-sw_init)
    append!(R, R_tmp)
    append!(pv_R, t_D_tmp)

    # saturation profile
    sor_ls = rel_perms_ls.sor
    sw_inj = core_flood.injected_water_saturation
    phi = core_props.porosity
    ut = core_flood.injection_velocity
    s1 = collect(linspace(min(sw_inj, 1-sor_ls-eps()), sw_shock_ls-eps(), 100))
    xt_s1 = dfw_ls.(s1)
    xt_shock_ls = 1/t_D_BT_ls
    xt_shock_hs = 1/t_D_BT_hs
    xt_prf=[xt_s1; xt_shock_ls; xt_shock_ls+eps(); xt_shock_hs; xt_shock_hs+eps(); 1/0.3]
    sw_prf=[s1; sw_shock_ls; sw_shock_hs; sw_shock_hs; 1-sor_hs; 1-sor_hs]

    return pv_R, R, xt_prf, sw_prf # for the time being to test the code
end

"""
the eq_const is defined as the ratio of the volume fraction of the solvent in oil to 
its volume fraction in water. For instance, For DME, the molar fraction is around 2.0
One criteria is that eq_const should probably be higher than 1.
"""
function water_soluble_solvent_flood(core_props, fluids_ls, fluids_hs, rel_perms_hs, 
    rel_perms_ls, core_flood, eq_const)
    # construct the fractional flow curves
    fw_ls, dfw_ls = fractional_flow_function(rel_perms_ls, fluids_ls)
    fw_hs, dfw_hs = fractional_flow_function(rel_perms_hs, fluids_hs)
    # low sal shock (tangent line from (0,0))
    point1 = (-eq_const/(1-eq_const), -eq_const/(1-eq_const))
    sw_shock_ls = tangent_line_saturation(rel_perms_ls, fluids_ls, point1)
    println("low sal sw_shock = $sw_shock_ls")
    t_D_BT_ls = 1/dfw_ls(sw_shock_ls) # breakthrough (BT) time [#PV]
    println("low sal breakthrough time = $t_D_BT_ls")
    # High sal shock (cross point between the ls tangent and the hs fw)
    sw_shock_hs = cross_point_saturation(fw_hs, rel_perms_hs, point1, (sw_shock_ls, fw_ls(sw_shock_ls)))
    println("high sal sw_shock = $sw_shock_hs")
    sw_init = core_flood.initial_water_saturation    
    t_D_BT_hs = (sw_shock_hs-sw_init)/(fw_hs(sw_shock_hs)-fw_hs(sw_init)) # breakthrough (BT) time [#PV]
    println("high sal breakthrough time = $t_D_BT_hs")
    
    # construct the recovery factor curve versus the # of PV
    R = zeros(1)
    pv_R = zeros(1)
    # at breakthrough of the hs brine
    push!(R, (1-fw_hs(sw_init))*t_D_BT_hs/(1-sw_init)) # recovery at BT
    push!(pv_R, t_D_BT_hs) # BT time
    # at breakthrough of the ls brine
    push!(R, R[end]+(1-fw_hs(sw_shock_hs))*(t_D_BT_ls-t_D_BT_hs)/(1-sw_init)) # recovery at BT
    push!(pv_R, t_D_BT_ls) # BT time

    # after breakthrough
    pv_inj = max(core_flood.injected_pore_volume, 2.0) # at least inject 2 pv
    f_sw = sw -> (pv_inj-1/dfw_ls(sw)) # find the outlet saturation at the end of injection
    sw_max = fzero(f_sw, sw_shock_ls)
    sw_tmp = linspace(sw_shock_ls, sw_max, 100)
    t_D_tmp = 1./dfw_ls.(sw_tmp)
    
    s_av_tmp = sw_tmp-(fw_ls.(sw_tmp)-1).*t_D_tmp
    R_tmp = (s_av_tmp-sw_init)/(1-sw_init)
    append!(R, R_tmp)
    append!(pv_R, t_D_tmp)

    # saturation profile
    sor = rel_perms_ls.sor
    sw_inj = core_flood.injected_water_saturation
    phi = core_props.porosity
    ut = core_flood.injection_velocity
    s1 = collect(linspace(min(sw_inj, 1-sor-eps()), sw_shock_ls-eps(), 100))
    xt_s1 = dfw_ls.(s1)
    xt_shock_ls = 1/t_D_BT_ls
    xt_shock_hs = 1/t_D_BT_hs
    xt_prf=[xt_s1; xt_shock_ls; xt_shock_ls+eps(); xt_shock_hs; xt_shock_hs+eps(); 1/0.3]
    sw_prf=[s1; sw_shock_ls; sw_shock_hs; sw_shock_hs; sw_init; sw_init]

    return pv_R, R, xt_prf, sw_prf # for the time being to test the code

end

end # module