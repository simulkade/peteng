# BL solution
# plotting the results of water-flooding oil recovery
# optimizing the rel-perm parametrs using the recovery data
using Roots, Dierckx, Plots
# include("rel_perms.jl")


# plot the fractional flow functions for a water-wet and an oil-wet system
include("rel_perms_real.jl")
# Water-wet Mixed-wet Oil-wet
# Swi 0.2 0.2 0.2
# Sor 0.1 0.2 0.3
# Ew=Krw(Sor) 0.6 0.8 0.8
# Eo=Kro(Swi) 0.8 0.8 0.5
# Nw 3 2 3
# No 2 2 3
krw0_ww = 0.3
krw0_ow = 1.0
kro0_ww = 0.6
kro0_ow = 0.76
nw_ww = 2.4
nw_ow= 2.4
no_ww = 2.0
no_ow= 2.0
sor_ww=0.1
sor_ow=0.12
swc_ww=0.09
swc_ow=0.09
SF=0.0 # 1.0 is water-wet, 0.0 is oil-wet
sw_plot_ow=linspace(swc_ow,1-sor_ow,100)
sw_plot_ww=linspace(swc_ww,1-sor_ww,100)

plot(sw_plot_ow, krw.(sw_plot_ow, krw0_ow, sor_ow, swc_ow, nw_ow),
  xlims=(0.0,1.0), linewidth=3, linestyle=:dash, linecolor=:blue,
  label="krw, water-wet")
plot!(sw_plot_ow, kro.(sw_plot_ow, kro0_ow, sor_ow, swc_ow, no_ow),
  xlims=(0.0,1.0), linewidth=3, linestyle=:solid, linecolor=:red,
  label="kro, water-wet")
plot!(sw_plot_ww, krw.(sw_plot_ww, krw0_ww, sor_ww, swc_ww, nw_ww),
  xlims=(0.0,1.0), linewidth=1, linestyle=:dash, linecolor=:blue,
  label="krw, oil-wet")
plot!(sw_plot_ww, kro.(sw_plot_ww, kro0_ww, sor_ww, swc_ww, no_ww),
  xlims=(0.0,1.0), linewidth=1, linestyle=:solid, linecolor=:red,
  label="kro, oil-wet")


# plot the recovery factors for a water-wet, a mixed-wet, and an oil-wet system
include("frac_flow_funcs.jl")
(xt_shock, sw_shock, xt_prf, sw_prf, t, p_inj, R_oil)=frac_flow_wf(muw=1e-3,
  muo=2e-3, ut=1e-4, phi=0.2,
  k=1e-12, swc=0.09, sor=0.12, kro0=0.76, no=2.0, krw0=1.0,
  nw=2.4, sw0=0.09, sw_inj=1.0, L=1, pv_inj=5)

# plot(t, R_oil)
