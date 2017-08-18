# plot residual oil saturation versus contact angle
using Plots
sor_ww=0.15
sor_ow=0.2
sor_mw=0.08
teta_ow=π
teta_ww=0.0
teta_mw=π/2+0.2 # just to shift it a bit
a=(sor_mw-0.5*(sor_ww-sor_ow)*cos(teta_mw)-0.5*(sor_ww+sor_ow))/(cos(teta_mw)^2-1)
b=0.5*(sor_ww-sor_ow)
c=0.5*(sor_ww+sor_ow)-a

teta=linspace(0, π, 100)
sor=a*cos(teta).^2+b*cos(teta)+c

plot(
plot(teta, sor, xlabel="contact angle [rad]", ylabel="Sor"),
plot(cos(teta), sor, xlabel="contact angle [rad]", ylabel="Sor")
)
