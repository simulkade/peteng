using Roots, Dierckx, NLopt, CoolProp, PyPlot
import JSON
include("rel_perms.jl")
include("frac_flow_funcs.jl")


# read the data files
lowsal=JSON.parsefile("../juliaphreeqc/data/Lowsal_data.json")
paper_title = "fathi2010"
exp_type = "coreFlood"
exp_name = ["SF20-a", "SF20-b"] # formation water and sea water
recovery_data = "coreFloodData" # -> time(day), recovery
core_flood_condition = lowsal["fathi2010"]["coreFlood"]["SF20-a"]
L_core = core_flood_condition["L(mm)"]/1000
D_core = core_flood_condition["D(mm)"]/1000
poros  = core_flood_condition["porosity"]/100
swi     = core_flood_condition["swi"]/100
soi     = core_flood_condition["soi"]/100 # 1-swi
T_c     = core_flood_condition["T(C)"]
oil_type = core_flood_condition["oil"]
mu_oil_25 = lowsal["fathi2010"]["oil"][oil_type]["viscosity(cp)"]/100 # Pa.s
t_day = lowsal["fathi2010"]["coreFloodData"]["SF20-a"]["time(day)"]
R_oil = lowsal["fathi2010"]["coreFloodData"]["SF20-a"]["recovery"]/100
injection_rate_pv_day = core_flood_condition["injectionRate(PV/day)"]
recovery_final = core_flood_condition["recovery(%OOIP)"]/100
dp_exp = 0.0

so_final = (1-recovery_final)*soi # final oil saturation
pv_inj = t_day[end]*injection_rate_pv_day
t_sec = t_day*24*3600 # s
# conversion of the input data to the right unit
A_core  = π*D_core^2/4
PV_core = L_core*A_core*poros
u_inj   = injection_rate_pv_day/(24*3600)*PV_core/A_core # Darcy velocity
T_K = T_c + 273.15
perm_ave = 0.002e-12 # m^2 permeability
mu_water = PropsSI("viscosity", "T", T_K, "P", 5e5, "H2O")
# My formulation for oil viscosity at higher Temperature
# source: P.Daučík et al., Temperature Dependence of the Viscosity of Hydrocarbon Fractions
# Acta Chimica Slovaca, Vol.1, No. 1, 2008, 43 – 57
# log(mu/mu0) = B(1/T - 1/T0)
# find B for n-Dodecane, and use it for extrapolation
mu_1   = PropsSI("viscosity", "T", 300, "P", 5e5, "n-Dodecane")
mu_2   = PropsSI("viscosity", "T", T_K, "P", 5e5, "n-Dodecane")
B_visc = log(mu_1/mu_2)/(1/300-1/T_K)
mu_oil = exp(B_visc*(1/T_K-1/(25+273.15)))*mu_oil_25 # Pa.s

# Corrections to make the data consistent;
# The recovery line must be a straight line at the begining of the experiment.
# I shift the time axis a bit to make it a straight line:
# plot(t_sec, R_oil, "-o")
(a,b) = linreg(t_sec[2:6], R_oil[2:6])
# 0 = a + b*t_shift
t_shift = -a/b
t_sec_cor = [0; t_sec[2:end]-t_shift]
# plot(t_sec_cor, R_oil, "-o")

# define the objective functions
# (xt_shock_ww, sw_shock_ww, xt_prf_ww, sw_prf_ww, t_ww, p_inj_ww, R_oil_ww)=
# frac_flow_wf(muw=mu_water, muo=mu_oil, ut=inj_vel, phi=poros_res,
#   k=perm_res, swc=swc_ww, sor=sor_ww, kro0=kro0_ww, no=no_ww, krw0=krw0_ww,
#   nw=nw_ww, sw0=swc_ww, sw_inj=1.0, L=L_res, pv_inj=inj_pv)


# define the objective functions
# param = [sor, swc, kro0, krw0, no, nw]
function error_calc(param, muw, muo, ut, phi, k, sw0, sw_inj, L, pv_inj, t_exp, R_exp, dp_exp=0.0)
  (sor, swc, kro0, krw0, no, nw) = param
  (xt_shock, sw_shock, xt_prf, sw_prf, t, p_inj, R)=
  frac_flow_wf(muw=muw, muo=muo, ut=ut, phi=phi,
    k=k, swc=swc, sor=sor, kro0=kro0, no=no, krw0=krw0,
    nw=nw, sw0=sw0, sw_inj=sw_inj, L=L, pv_inj=pv_inj)
  R_int = Spline1D(t, R, k=1)
  dp_int = Spline1D(t, p_inj, k=1)
  error_vals_R = R_int(t_exp) - R_exp
  error_vals = error_vals_R
  if dp_exp != 0.0
      error_vals_dp = (dp_int(t_exp) - dp_exp)/maximum(dp_exp) # normalize
      error_vals = [error_vals_R; error_vals_dp;]
  end
  return error_vals
end

function error_calculation(param, param_ind, muw, muo, ut, phi, k, sw0, sw_inj, L, pv_inj,
  t_exp, R_exp, dp_exp; grad_calc=true)
  eps1 = 1e-8
  error_vals = error_calc(param, muw, muo, ut, phi, k, sw0, sw_inj, L, pv_inj, t_exp, R_exp, dp_exp)
  if grad_calc
    jac  = zeros(length(error_vals), length(param_ind))
  else
    jac = zeros(0)
  end
  if grad_calc
    for j in eachindex(param_ind)
      param2 = copy(param)
      param2[param_ind[j]]+=eps1
      error_val2 = error_calc(param2, muw, muo, ut, phi, k, sw0, sw_inj, L, pv_inj, t_exp, R_exp, dp_exp)
      jac[:, j] = (error_val2 - error_vals) / eps1
    end
  end
  return error_vals, jac
end


function objective_function(param::Vector{Float64}, grad::Vector{Float64}, param_all,
  param_ind, muw, muo, ut, phi, k, sw0, sw_inj, L, pv_inj, t_exp, R_exp, dp_exp; w=1.0)
  param2 = copy(param_all)
  param2[param_ind] = param
  grad_calc = false
  if length(grad)>0
    grad_calc = true
  end
  error_val, jac = error_calculation(param2, param_ind, muw, muo, ut, phi, k, sw0,
    sw_inj, L, pv_inj, t_exp, R_exp, dp_exp, grad_calc = grad_calc)
  obj_func_value = sum(abs2, w.*error_val)
  if length(grad)>0
      grad[:] = 2.0*sum(w.*error_val.*jac, 1)
  end
  return obj_func_value
end

function plot_results(param, muw, muo, ut, phi, k, sw0, sw_inj, L, pv_inj, t_exp, R_exp, dp_exp)
  (sor, swc, kro0, krw0, no, nw) = param
  (xt_shock, sw_shock, xt_prf, sw_prf, t, p_inj, R)=
  frac_flow_wf(muw=muw, muo=muo, ut=ut, phi=phi,
    k=k, swc=swc, sor=sor, kro0=kro0, no=no, krw0=krw0,
    nw=nw, sw0=sw0, sw_inj=sw_inj, L=L, pv_inj=pv_inj)
  R_int = Spline1D(t, R, k=1)
  error_vals = R_int(t_exp) - R_exp
  plot(t_exp, R_int(t_exp), t_exp, R_exp, "o")
end

plot_quick = param -> plot_results(param, mu_water, mu_oil, u_inj, poros, perm_ave, swi, 1.0, L_core, pv_inj, t_sec_cor, R_oil, dp_exp)
# Optimization problem:

# parameter values
# [sor, swc, kro0, krw0, no, nw]
param_all = [so_final, 0.2, 0.5, 0.5, 2, 2.4]
param_ind = [1, 2,3,4,5,6]
w = ones(length(R_oil))
w[end]=2
w[end-1]=2
w[7] = 3
obj_fun = (param, grad) -> objective_function(param::Vector{Float64}, grad::Vector{Float64}, param_all,
  param_ind, mu_water, mu_oil, u_inj, poros, perm_ave, swi, 1.0, L_core, pv_inj, t_sec_cor, R_oil, dp_exp, w=w)

# x_init = [0.2, 0.8, 0.4, 2, 2]
# x_lb = [0.1, 0.1, 0.1, 0.5, 0.5]
# x_ub = [0.5, 1.0, 1.0, 5.0, 5.0]

x_init = copy(param_all)
x_lb = [0.08, 0.08, 0.05, 0.05, 1.0, 1.0]
x_ub = [0.5, 0.5, 1.0, 1.0, 4.0, 4.0]

# algorithms
# :LD_MMA
# :LN_COBYLA
# :LD_LBFGS
# :GN_DIRECT
# :GN_DIRECT_L
opt_alg=:LD_MMA

opt1 = Opt(opt_alg, length(x_init)) # choose the algorithm
lower_bounds!(opt1, x_lb)
upper_bounds!(opt1, x_ub)
ftol_rel!(opt1, 1e-10)
ftol_abs!(opt1, 1e-10)

min_objective!(opt1, obj_fun)
(fObjOpt, paramOpt, flag) = optimize(opt1, x_init)

param_plot = copy(param_all)
param_plot[param_ind] = paramOpt
plot_quick(param_plot)
