# BL solution
# plotting the results of water-flooding oil recovery
# optimizing the rel-perm parametrs using the recovery data
using Roots, Dierckx, Plots, JLD
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
sor_ww=0.15
sor_ow=0.2
sor_mw=0.08
swc_ww=0.09
swc_ow=0.09
teta_ow=deg2rad(180)
teta_ww=deg2rad(0)
teta_mw=deg2rad(90)
teta_prime=teta_mw-π/2 # see my report
a=(sor_mw-0.5*(sor_ww-sor_ow)*cos(teta_mw)-0.5*(sor_ww+sor_ow))/(cos(teta_mw)^2-1)
b=0.5*(sor_ww-sor_ow)
c=0.5*(sor_ww+sor_ow)-a

SF=0.0 # 1.0 is water-wet, 0.0 is oil-wet
sw_plot_ow=linspace(swc_ow,1-sor_ow,100)
sw_plot_ww=linspace(swc_ww,1-sor_ww,100)

plot(size=(500,400), xtickfont = font(12, "Courier"), ytickfont=font(12, "Courier"),
  ylabel="krw, kro", xlabel="Water saturation", legendfont=font(12, "Courier"),
  guidefont=font(12, "Courier"))
plot!(sw_plot_ow, krw.(sw_plot_ow, krw0_ow, sor_ow, swc_ow, nw_ow),
  xlims=(0.0,1.0), linewidth=3, linestyle=:dash, linecolor=:blue,
  label="krw, Oil-wet")
plot!(sw_plot_ow, kro.(sw_plot_ow, kro0_ow, sor_ow, swc_ow, no_ow),
  xlims=(0.0,1.0), linewidth=3, linestyle=:solid, linecolor=:red,
  label="kro, Oil-wet")
plot!(sw_plot_ww, krw.(sw_plot_ww, krw0_ww, sor_ww, swc_ww, nw_ww),
  xlims=(0.0,1.0), linewidth=1, linestyle=:dash, linecolor=:blue,
  label="krw, Water-wet")
plot!(sw_plot_ww, kro.(sw_plot_ww, kro0_ww, sor_ww, swc_ww, no_ww),
  xlims=(0.0,1.0), linewidth=1, linestyle=:solid, linecolor=:red,
  label="kro, Water-wet")
savefig("two_rel_perms.png")

# plot the recovery factors for a water-wet, a mixed-wet, and an oil-wet system
include("frac_flow_funcs.jl")
L_res=1000 # m
perm_res=0.001e-12 # Darcy
inj_pv=5.0 # pore volumes
poros_res=0.35 # reservoir porosity
inj_vel=1.0/(3600*24) # m/s
mu_oil=2e-3 # oil viscosity
mu_water=1e-3 # water viscosity

# water-wet reservoir
(xt_shock_ww, sw_shock_ww, xt_prf_ww, sw_prf_ww, t_ww, p_inj_ww, R_oil_ww)=
frac_flow_wf(muw=mu_water, muo=mu_oil, ut=inj_vel, phi=poros_res,
  k=perm_res, swc=swc_ww, sor=sor_ww, kro0=kro0_ww, no=no_ww, krw0=krw0_ww,
  nw=nw_ww, sw0=swc_ww, sw_inj=1.0, L=L_res, pv_inj=inj_pv)
# oil-wet reservoir
(xt_shock_ow, sw_shock_ow, xt_prf_ow, sw_prf_ow, t_ow, p_inj_ow, R_oil_ow)=
frac_flow_wf(muw=mu_water, muo=mu_oil, ut=inj_vel, phi=poros_res,
  k=perm_res, swc=swc_ow, sor=sor_ow, kro0=kro0_ow, no=no_ow, krw0=krw0_ow,
  nw=nw_ow, sw0=swc_ow, sw_inj=1.0, L=L_res, pv_inj=inj_pv)

plot(t_ow/(24*3600), R_oil_ow, size=(400,300), linecolor=:black, linewidth=2,
  xtickfont = font(10, "Courier"), ytickfont=font(10, "Courier"),
  xlabel="time [day]", ylabel="Recovery factor", label="Oil-wet",
  legendfont=font(10, "Courier"), guidefont=font(12, "Courier"), ylims=(0.0,1.0))
# savefig("recovery_time_ow.png")
plot!(t_ww/(24*3600), R_oil_ww, linecolor=:blue, label="Water-wet", linewidth=2,
  linestyle=:dash)
# savefig("recovery_time_ww_ow.png")

# BL for different wettability
plot(size=(500,400), xtickfont = font(10, "Courier"), ytickfont=font(10, "Courier"),
  xlabel="t [day]", ylabel="Recovery factor", legendfont=font(10, "Courier"),
  guidefont=font(12, "Courier"))
for SF in 0.0:0.25:1.0
  krw0=krw0_ww*SF+krw0_ow*(1-SF)
  kro0=kro0_ww*SF+kro0_ow*(1-SF)
  #sor=sor_ww*SF+sor_ow*(1-SF)
  swc=swc_ww*SF+swc_ow*(1-SF)
  no= no_ww*SF+no_ow*(1-SF)
  nw= nw_ww*SF+nw_ow*(1-SF)
  teta=teta_ww*SF+teta_ow*(1-SF)
  # sor=(max(sor_ow, sor_ww)-sor_mw)*cos(teta-teta_prime)^2+sor_mw
  sor=a*cos(teta).^2+b*cos(teta)+c
  (xt_shock, sw_shock, xt_prf, sw_prf, t, p_inj, R_oil)=
  frac_flow_wf(muw=mu_water, muo=mu_oil, ut=inj_vel, phi=poros_res,
    k=perm_res, swc=swc, sor=sor, kro0=kro0, no=no, krw0=krw0,
    nw=nw, sw0=swc, sw_inj=1.0, L=L_res, pv_inj=inj_pv)
    plot!(t/(24*3600), R_oil, size=(500,400), linewidth=2,
      xtickfont = font(10, "Courier"), ytickfont=font(10, "Courier"),
      xlabel="time [day]", ylabel="Recovery factor", label="$SF Sal. Factor",
      legendfont=font(10, "Courier"), guidefont=font(12, "Courier"), ylims=(0.0,1.0))
    savefig("recovery_forced_imb_$SF-WW.png")
end
plot!()


# plot(t, R_oil)
# plot the imbibition Pc-Sw curve for different SF
plot(size=(500,400), xtickfont = font(10, "Courier"), ytickfont=font(10, "Courier"),
  xlabel="Sw", ylabel="Pc [Pa]", legendfont=font(10, "Courier"),
  guidefont=font(12, "Courier"))
for SF in 0.1:0.2:0.9
  krw0=krw0_ww*SF+krw0_ow*(1-SF)
  kro0=kro0_ww*SF+kro0_ow*(1-SF)
  #sor=sor_ww*SF+sor_ow*(1-SF)
  swc=swc_ww*SF+swc_ow*(1-SF)
  no= no_ww*SF+no_ow*(1-SF)
  nw= nw_ww*SF+nw_ow*(1-SF)
  teta=teta_ww*SF+teta_ow*(1-SF)
  # sor=(max(sor_ow, sor_ww)-sor_mw)*cos(teta-teta_prime)^2+sor_mw
  sor=a*cos(teta).^2+b*cos(teta)+c
  pce= 5e2
  sw_plot=linspace(0.0,1.0,100)

  pc=pc_imb.(sw_plot, pce, swc, sor, teta=teta)
  plot!(sw_plot, pc, xlims=(0.0,1.0), ylims=(-2000, 2000), linewidth=2,
  label="$SF Sal. Factor")
  savefig("pc_sw_$SF-WW.png")
end
plot!()
