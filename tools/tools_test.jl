# test the tools module
include("tools.jl")
LI = LinearInterp

using Plots

x = collect(sort(rand(10)))
y = sin.(10*x)

xfine = collect(x[1]:0.01:x[end])

intp1 = LI.interp1D(x, y)

yfine = intp1(xfine)

scatter(x, y, marker=:circle)
scatter!(xfine, yfine, marker=:s)
