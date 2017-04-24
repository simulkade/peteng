# all the plots for my presentation in IOR2017, Stavanger
using Plots
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
krw0_ww = 0.33
krw0_ow = 0.28
kro0_ww = 0.79
kro0_ow = 0.72
nw_ww = 2.4
nw_ow= 2.4
no_ww = 2.0
no_ow= 1.6
sor_ww=0.145
sor_ow=0.37
sor_mw=0.08
swc_ww=0.08
swc_ow=0.08
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
  xlims=(0.0,1.0), linewidth=3, linestyle=:dash, linecolor=:red,
  label="krw, Oil-wet")
plot!(sw_plot_ow, kro.(sw_plot_ow, kro0_ow, sor_ow, swc_ow, no_ow),
  xlims=(0.0,1.0), linewidth=3, linestyle=:dash, linecolor=:red,
  label="kro, Oil-wet")
plot!(sw_plot_ww, krw.(sw_plot_ww, krw0_ww, sor_ww, swc_ww, nw_ww),
  xlims=(0.0,1.0), linewidth=1, linestyle=:solid, linecolor=:blue,
  label="krw, Water-wet")
plot!(sw_plot_ww, kro.(sw_plot_ww, kro0_ww, sor_ww, swc_ww, no_ww),
  xlims=(0.0,1.0), linewidth=1, linestyle=:solid, linecolor=:blue,
  label="kro, Water-wet")
savefig("rel_perms_fathi.png")
