module FractionalFlow

using Roots, Dierckx, PyPlot, JFVM, ProgressMeter, JSON
import Optim
# import GR
include("CoreyFunctions.jl")
include("water_flood.jl")
include("low_sal_flood.jl")
include("water_solvent_flood.jl")
include("tertiary_flood.jl")
include("water_flood_fvm.jl")
include("water_flood_fvm_upwind.jl")
include("water_flood_upwind_BL.jl")
include("fractional_flow_cases.jl")
include("HTMLElements.jl")
include("fractional_flow_input.jl")

CF = CoreyFunctions # a shorter name

# types:
"""
Core relative permeability parameters
    krw0::Real
    kro0::Real
    swc::Real
    sor::Real
    nw::Real
    no::Real
"""
struct CoreyRelativePermeability
    krw0::Real
    kro0::Real
    swc::Real
    sor::Real
    nw::Real
    no::Real
end

"""
The properties of the core are stored in this structure
    length::Real
    diameter::Real
    porosity::Real
    permeability::Real
"""
struct CoreProperties
    length::Real
    diameter::Real
    porosity::Real
    permeability::Real
    surface_area::Real
    matrix_density::Real
end

"""
species concentration
    species::AbstractString
    molality::Real
"""
struct SpeciesConcentration
    species::AbstractString
    concentration::Real
end

"""
Defines a brine solution, phreeqc style
    name::AbstractString
    number::Integer
    temperature::Real
    pressure::Real
    unit::AbstractString
    salinity::Array{SpeciesConcentration, 1}
    pH::Real
"""
struct Brine
    name::AbstractString
    number::Integer
    temperature::Real
    pressure::Real
    unit::AbstractString
    salinity::Array{SpeciesConcentration, 1}
    pH::Real
    density::Real
    viscosity::Real
end

struct Oil
    density::Real
    viscosity::Real
    acid_number::Real
    base_number::Real
    asphaltene::Real
end


"""
A structure for storing the core flooding experimental conditions
    injection_velocity::Real        # m/s
    injected_pore_volume::Real
    back_pressure::Real             # Pa
    initial_water_saturation::Real
    injected_water_saturation::Real
"""
struct CoreFlooding
    injection_velocity::Real # m/s
    injected_pore_volume::Real
    back_pressure::Real
    initial_water_saturation::Real
    injected_water_saturation::Real
end

"""
a structure that stores the viscosity of the fluids
"""
struct Fluids
    oil_viscosity::Real
    water_viscosity::Real
end

"""
A line to be used for storing the shock line results
"""
struct Line
    start_point::Array{Real,1}
    end_point::Array{Real,1}
end

"""
A structure that stores the fractional flow curves and the final solution
    fractional_flow_functions::Array{Function, 1}
    shock_line::Array{Real, 2}
    recovery_pv::Array{Real, 2}
    recovery_time::Array{Real, 2}
    saturation_profile_xt::Array{Real, 2}
"""
struct FracFlowResults
    fractional_flow_functions::Array{Function, 1}
    shock_lines::Array{Line, 1}
    recovery_pv::Array{Real, 2}
    recovery_time::Array{Real, 2}
    saturation_profile_xt::Array{Real, 2}
    tracer_profile_xt::Array{Real, 2}
    dp_pv::Array{Real, 2}
    dp_time::Array{Real, 2}
    water_cut_pv::Array{Real, 2}
    water_cut_time::Array{Real, 2}
end

# functions

"""
from old julia
"""
function linspace(a,b,c)
    return range(a, stop=b, length=c)
end

function indmin(a)
    return argmin(a)
end

function indmax(a)
    return argmax(a)
end


"""
visualize(rel_perm::CoreyRelativePermeability)
This function visualizes the relative permeability curve that is defined
by a CoreyRelativePermeability structure
"""
function visualize(rel_perms::CoreyRelativePermeability; label="")
    # convert the rel_perm structure to the rel_perm functions
    krw, kro, dkrwdsw, dkrodsw = rel_perm_functions(rel_perms)
    sw = linspace(0,1, 100)
    plot(sw, krw.(sw), label="krw-$label")
    plot(sw, kro.(sw), label="kro-$label")
    xlabel("Sw")
    ylabel("Relative permeability")
    legend()
end

"""
tabulate the fluids viscosity
"""
function print_fluids(fluids::Fluids; title = nothing)
    HTMLElements.table([fluids.oil_viscosity fluids.water_viscosity],
        header_row=[:mu_water, :mu_oil], title=title)
end

"""
tabulate the core properties
"""
function print_core_properties(core_props::CoreProperties; title=nothing)
    core_param = [core_props.length core_props.diameter core_props.porosity core_props.permeability]
    HTMLElements.table(core_param, header_row=[:length, :diameter, :porosity, :permeability], title=title)
end

"""
tabulate the ralative permeability table in a readable format
"""
function print_relperm(rel_perm::CoreyRelativePermeability; title=nothing)
    rel_perm_params = [rel_perm.krw0 rel_perm.kro0 rel_perm.nw rel_perm.no rel_perm.swc rel_perm.sor]
    HTMLElements.table(rel_perm_params, header_row = [:krw0, :kro0, :nw, :no, :Swc, :Sor], title=title)
end

"""
tabulates the core flooding conditions
"""
function print_core_flood(core_flood::CoreFlooding; title=nothing)
    core_flood_param = [core_flood.injection_velocity core_flood.injected_pore_volume core_flood.initial_water_saturation core_flood.injected_water_saturation]
    header_row = [:u_inj_m_s, :pv_injected, :Sw_init, :Sw_inj]
    HTMLElements.table(core_flood_param, header_row = header_row, title=title)
end


"""
This function visualizes the fractional flow curves
"""
function visualize(rel_perms::CoreyRelativePermeability, fluids::Fluids; label="")
    # convert the rel_perm structure to the rel_perm functions
    fw, dfw = fractional_flow_function(rel_perms, fluids)
    sw = linspace(0,1, 100)
    subplot(2,1,1)
    plot(sw, fw.(sw), label="fw-$label")
    xlabel("Sw [-]")
    ylabel("Fractional flow [-]")
    legend()
    subplot(2,1,2)
    plot(sw, dfw.(sw), label="dfw/dSw-$label")
    xlabel("Sw [-]")
    ylabel("Fractional flow derivative [-]")
    legend()
end

"""
Only visualizes the solution procedure
"""
function visualize_solution(wf_res::FracFlowResults)
    # plot the fractional flow functions and the solution procedure
    fw = wf_res.fractional_flow_functions
    sw = linspace(0, 1, 100)
    for f in fw
        plot(sw, f.(sw))
    end
    shock_lines = wf_res.shock_lines
    for sl in shock_lines
        plot(sl.start_point[1], sl.start_point[2], "ro")
        plot(sl.end_point[1], sl.end_point[2], "ro")
        plot([sl.start_point[1], sl.end_point[1]], [sl.start_point[2], sl.end_point[2]])
    end
    xlabel("Water saturation [-]")
    ylabel("water fractional flow [-]")
    legend()
end

"""
visualize the saturation and concentration profiles
"""
function visualize_profiles(wf_res::FracFlowResults)
    plot(wf_res.saturation_profile_xt[:,1], wf_res.saturation_profile_xt[:,2], label = "Sw")
    plot(wf_res.tracer_profile_xt[:,1], wf_res.tracer_profile_xt[:,2], label = "tracer")
    xlabel("xD/tD [-]")
    ylabel("Water saturation")
    legend()
end


"""
visualize the fractional flow results
"""
function visualize(wf_res::FracFlowResults)
    # plot the fractional flow functions and the solution procedure
    fw = wf_res.fractional_flow_functions
    sw = linspace(0, 1, 100)
    figure()
    for f in fw
        plot(sw, f.(sw))
    end
    shock_lines = wf_res.shock_lines
    for sl in shock_lines
        plot(sl.start_point[1], sl.start_point[2], "ro")
        plot(sl.end_point[1], sl.end_point[2], "ro")
        plot([sl.start_point[1], sl.end_point[1]], [sl.start_point[2], sl.end_point[2]])
    end
    xlabel("Water saturation [-]")
    ylabel("water fractional flow [-]")
    # plot the recovery factors (pv)
    figure()
    plot(wf_res.recovery_pv[:,1], wf_res.recovery_pv[:,2])
    xlabel("Pore volume [-]")
    ylabel("Recovery factor")
    figure()
    plot(wf_res.dp_pv[:,1], wf_res.dp_pv[:,2])
    xlabel("Pore volume [-]")
    ylabel("Pressure drop [Pa]")
    figure()
    plot(wf_res.water_cut_pv[:,1], wf_res.water_cut_pv[:,2])
    xlabel("Pore volume [-]")
    ylabel("Water cut [-]")
    # plot the recovery factors (time)
    # figure()
    # plot(wf_res.recovery_time[:,1], wf_res.recovery_time[:,2])
    # xlabel("time [s]")
    # ylabel("Recovery factor")
    # plot the saturation profile
    figure()
    plot(wf_res.saturation_profile_xt[:,1], wf_res.saturation_profile_xt[:,2], label = "Sw")
    plot(wf_res.tracer_profile_xt[:,1], wf_res.tracer_profile_xt[:,2], label = "tracer")
    xlabel("xD/tD [-]")
    ylabel("Water saturation")
    legend()
end

"""
corey_rel_perm = oil_water_rel_perms(krw0=0.4, kro0=0.9,
    swc=0.15, sor=0.2, nw=2.0, no = 2.0)
This functions defines a CoreyRelativePermeability structure
with a set of predefined values
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

function core_properties(;L=0.15, D=0.03, φ=0.3, k=1e-12, a=2000, ρ=2700)
    return CoreProperties(L, D, φ, k, a, ρ)
end

"""
rel_perm_functions(rel_perm::CoreyRelativePermeability)
returns:
     krw, kro, dkrwdsw, dkrodsw
"""
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
    eps1 = 1e-5 # distance from the end point saturations
    swc = rel_perms.swc
    sor = rel_perms.sor
    sw0 = sw_fw_point[1]
    fw0 = sw_fw_point[2]
    fw, dfw = fractional_flow_function(rel_perms, fluids)

    f_shock = sw -> (1-((fw(sw)-fw0)/(sw-sw0))/dfw(sw))
    sw_tmp = collect(linspace(swc+eps1, 1-sor-eps1, 1000))
    # ind_min = indmin(abs.(f_shock.(sw_tmp)))
    ind_max = indmax(dfw.(sw_tmp))
    sw_max = sw_tmp[ind_max]
    res = Optim.optimize(x->-dfw(x), swc, 1-sor, Optim.GoldenSection())
    sw_max = res.minimizer
    sw_tmp = linspace(sw_max, 1-sor-eps1, 100) # sw_tmp[ind_max:end]
    f_tmp = f_shock.(sw_tmp)
    # println(f_tmp)
    sw_right = 1-sor-eps1
    try
        sw_right = find(f_tmp.<0.0)[1]
    catch
        sw_right = 1-sor-eps1
        # @info "difficulty finding the shock front saturation range!"
    end
    ind_min = indmin(abs.(f_shock.(sw_tmp)))
    # ind_shock = indmax(abs.((fw.(sw_tmp[ind_max:end])-fw0)./(sw_tmp[ind_max:end]-sw0)))
    sw_shock_est = sw_tmp[ind_min]
    # println(sw_shock_est)
    sw_shock = 0.0
    try
        if f_shock(sw_max)*f_shock(sw_right)>0
            sw_shock = fzero(f_shock, sw_shock_est) #  ftol = 1e-10, xtol = 1e-7)
            if (abs(f_shock(sw_shock))>abs(f_shock(sw_shock_est))) || (sw_shock<sw_max) || (sw_shock>1-sor-eps1)
                sw_shock = sw_shock_est
            end
        else
            sw_shock = fzero(f_shock, [sw_max,sw_right])
            if abs(f_shock(sw_shock))>abs(f_shock(sw_shock_est))
                sw_shock = sw_shock_est
            end
        end
    catch
        sw_shock = sw_shock_est
        # @info "shock front saturation is estimated: $sw_shock, error is $(f_shock(sw_shock))"
    end
    # println(abs(f_shock(sw_shock_est)))
    # println(abs(f_shock(sw_shock)))
    return sw_shock
end

function cross_point_saturation(fw, rel_perms, point1, point2; sw_left=0, sw_right=0)
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
    if sw_right==0
        sw_left, sw_right = minmax(sw_neg, sw_pos)
    end
    sw_cross = sw_tmp[indmin(abs.(f_tmp))]
    try
        if f(sw_left)*f(sw_right)>0
            sw_cross = fzero(f, sw_cross)
        else
            sw_cross = fzero(f, [sw_left, sw_right])
        end
    catch
        # @info "Cross point is estimated: $sw_cross, error is $(f(sw_cross))"
    end
    return sw_cross
end

"""
fint the saturation at the outlet of the core at the end of injection
"""
function outlet_saturation(pv_inj, sw_shock, dfw, rel_perms)
    eps1 = 1e-2
    sor = rel_perms.sor
    sw_tmp = collect(linspace(sw_shock+eps(), 1-sor-eps1, 1000))
    f = sw -> (pv_inj-1/dfw(sw)) # find the outlet saturation at the end of injection
    f_tmp = f.(sw_tmp)
    sw_left = sw_shock+eps()
    sw_right = 1-sor-eps1
    # println(f_tmp)
    sw_max = 0.0
    try
        sw_max = sw_tmp[find(f_tmp.<0.0)[1]]
    catch
        sw_max = sw_tmp[indmin(f_tmp)]
    end
    try
        if f(sw_left)*f(sw_right)>0
            sw_max = fzero(f_sw, [sw_left, sw_right])
        else
            sw_max = fzero(f, [sw_left, sw_right])
        end
    catch
        sw_max = sw_tmp[indmin(abs.(f_tmp))]
        # @info "Outlet saturation is estimated: $sw_max, error is $(f(sw_max))"
    end
    return sw_max
end

function trapz(x,y)
    # copied from https://github.com/hwborchers/NumericalMath.jl/blob/master/src/integrate.jl
        n=length(x)
        if (length(y) != n)
            error("Vectors 'x', 'y' must be of same length")
        end
        r=0
        for i in 2:n
            r+=(x[i] - x[i-1])*(y[i] + y[i-1])
        end
        return r/2.0
    end

end # module
