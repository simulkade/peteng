# load and plot the imbibition data
include("rel_perms.jl")
include("imbibition_IMPES.jl")

plot(size=(500,400), xtickfont = font(10, "Courier"), ytickfont=font(10, "Courier"),
  xlabel="t [day]", ylabel="Recovery Factor", legendfont=font(10, "Courier"),
  guidefont=font(12, "Courier"))
for SF in 0.1:0.2:0.9
  (t,R)=imb_impes(SF)
  save("imb_$SF-ww.jld", "t", t, "R", R)
  plot!(t/(24*3600), R, xlims=(0.0,50.0), ylims=(0, 1), linewidth=2,
  label="$SF water-wet")
  savefig("imb_recovery_$SF-WW.png")
end
plot!()
