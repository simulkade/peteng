using Roots, Dierckx, NLopt, CoolProp, PyPlot, JLD, PyCall
import JSON
@pyimport scipy.optimize as sciopt
include("../functions/rel_perms_real.jl")
include("../functions/frac_flow_funcs.jl")


# read the data files
L_core = 0.05
D_core = 0.03
poros  = 0.256
A_core  = π*D_core^2/4
PV_core = L_core*A_core*poros

swi     = 0.334
soi     = 1-swi
T_c     = 55.6
mu_oil = 0.35e-3
pv_inj_R_exp = [0.0, 0.029, 0.152, 0.274, 1.195, 2.177, 3.527, 4.939, 6.289, 7.700, 9.051]
R_oil_exp = [0.0, 0.043, 0.236, 0.430, 0.534, 0.564, 0.584, 0.593, 0.607, 0.610, 0.610]
R_oil_exp_int = Spline1D(pv_inj_R_exp, R_oil_exp, k=1)
pv_inj_dp = [0.08, 0.18, 0.24, 0.43, 0.61, 0.86, 1.43, 2.20, 3.02, 3.96, 4.90,
5.84, 6.80, 7.73, 8.73]
dp_exp = [0.001098, 0.486813, 0.53956, 0.486813, 0.434065, 0.381318, 0.334065, 0.302197,
0.271428, 0.261538, 0.247252, 0.23956, 0.231868, 0.224175, 0.221978]*20*1e5

pv_inj_exp = pv_inj_dp[:] # same as the time steps for the dp
R_oil = R_oil_exp_int(pv_inj_exp) # just to have both dp and R for the same time steps

pv_inj = pv_inj_exp[end]
injection_rate = 0.001/(3600*24) # m^3/day
recovery_final = R_oil_exp[end]

so_final = (1-recovery_final)*soi # final oil saturation

t_sec = pv_inj_exp*PV_core/injection_rate

# conversion of the input data to the right unit
u_inj   = injection_rate/A_core # Darcy velocity
T_K = T_c + 273.15
perm_ave = 0.001e-12 # m^2 permeability
mu_water = PropsSI("viscosity", "T", T_K, "P", 5e5, "H2O")

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
  dp_int = Spline1D(t, p_inj, k=1)
  error_vals = R_int(t_exp) - R_exp
  t_plot = linspace(0.0, maximum(t_exp), 100)
  figure(1)
  plot(t_plot, R_int(t_plot), t_exp, R_exp, "o")
  xlabel("time [s]")
  ylabel("Recovery factor")
  figure(2)
  plot(t_plot, dp_int(t_plot))
  if dp_exp !=0.0
      plot(t_exp, dp_exp, "o")
  end
  xlabel("time [s]")
  ylabel("pressure drop [Pa]")
end

plot_quick = param -> plot_results(param, mu_water, mu_oil, u_inj, poros, perm_ave, swi, 1.0, L_core, pv_inj, t_sec, R_oil, dp_exp)
# Optimization problem:

# parameter values
# [sor, swc, kro0, krw0, no, nw]
param_all = [0.2, 0.3, 1.0, 1.0, 3, 3]
param_ind = [1, 2,3,4,5,6]
w = ones(2*length(R_oil))
# w[end]=2
# w[end-1]=2
# w[7] = 3
obj_fun = (param, grad) -> objective_function(param::Vector{Float64}, grad::Vector{Float64}, param_all,
  param_ind, mu_water, mu_oil, u_inj, poros, perm_ave, swi, 1.0, L_core, pv_inj, t_sec, R_oil, dp_exp, w=w)

# x_init = [0.2, 0.8, 0.4, 2, 2]
# x_lb = [0.1, 0.1, 0.1, 0.5, 0.5]
# x_ub = [0.5, 1.0, 1.0, 5.0, 5.0]

x_init = param_all[:]
x_lb = [0.05, 0.25, 0.05, 0.05, 1.0, 1.0]
x_ub = [0.45,0.45, 1.0, 1.0, 4.0, 4.0]

# algorithms
# :LD_MMA
# :LN_COBYLA
# :LD_LBFGS
# :GN_DIRECT
# :GN_DIRECT_L
opt_alg=:LN_COBYLA

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

# optimized parameters:
# 0.24774  0.315125  0.948653  0.972952  2.60396  2.92036
