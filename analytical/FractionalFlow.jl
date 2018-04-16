module FractionalFlow

using Roots, Dierckx, PyPlot
include("CoreyFunctions.jl")
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
    if !isempty(rel_perms)
        if rel_perms.swc>sw_init
            swc = rel_perms.swc
        end
    end
    return CoreFlooding(u_inj, pv_inject, p_back, swc, sw_inj)
end

function core_properties(;L=0.15, D=0.03, φ=0.3, k=1e-12)
    return CoreProperties(L, D, φ, k)
end

function test_water_flood()
    fluids = oil_water_fluids()
    rel_perms = oil_water_rel_perms()
    core_flood = core_flooding()
    core_props = core_properties()
    pv, R, xt, sw = water_flood(core_props, fluids, rel_perms, core_flood)
    fw, dfw = fractional_flow_function(rel_perms, fluids)
    sw_tmp = linspace(0,1,100)
    # plot(sw_tmp, fw.(sw_tmp), xlabel = "Sw", ylabel="fw", label="")
    # plot!(sw_tmp, dfw.(sw_tmp))
    figure(1)
    plot(pv, R) 
    xlabel("PV injected")
    ylabel("Recovery factor") 
    title("Water flooding")
    figure(2)
    plot(xt, sw)
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
    m = (point2[2]-point1[2])/(point2[1]-point1[1])
    b = point1[2]-m*point1[1]
    f = sw -> m*sw+b-fw(sw)
    sw_cross = fzero(f, [swc+eps(), 1-sor-eps()])
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
    s1 = collect(linspace(min(sw_inj, 1-sor-eps()), sw_shock, 100))
    xt_s1 = ut/phi*dfw.(s1)
    xt_shock = ut/phi*dfw(sw_shock)
    xt_prf=[xt_s1; xt_shock; xt_shock+eps(); 2*xt_shock]
    sw_prf=[s1; sw_shock; sw_init; sw_init]

    return pv_R, R, xt_prf, sw_prf # for the time being to test the code

end

function low_sal_water_flood(core_props, fluids, rel_perms_hs, 
    rel_perms_ls, core_flood)

end

function single_ion_adsorption_water_flood()

end

function water_soluble_solvent_flood()

end

end # module