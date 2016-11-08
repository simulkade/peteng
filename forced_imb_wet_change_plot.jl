# imbibition with wettability change
# solve transport for one species in the aqueous phase
# change the SF accordingly
# plot the production curve
import Plots
include("forced_imbibition_IMPES.jl")
sw0=0.1
sw=0.1
Plots.plot(size=(500,400), xtickfont = Plots.font(10, "Courier"), ytickfont=Plots.font(10, "Courier"),
  xlabel="t [s]", ylabel="Recovery factor", legendfont=Plots.font(10, "Courier"),
  guidefont=Plots.font(12, "Courier"))
t_plot=0.0
R_plot=0.0
for SF0 in 0.1:0.2:0.9
  if sw!=sw0
    sw0=sw.value[2:end-1]
  end
  t,R,sw=forced_imb_impes(SF0, sw0)
  t_plot=t_plot[end]+t
  R_plot=R_plot[end]+R*(1-R_plot[end])
  Plots.plot!(t_plot, R_plot, ylims=(0.0,1.0), linewidth=2,
  label="$SF0 Sal. Factor")
end
Plots.plot!()
Plots.savefig("forced_imb_multiple_sal.png")
