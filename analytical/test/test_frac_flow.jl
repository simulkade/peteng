# testing some of the functionalities of the fractional flow package

include("../FractionalFlow/FractionalFlow.jl")
FF = FractionalFlow

## 1- read the input file
core_props, fluids, rel_perms, core_flood = FF.read_frac_flow_input("input_file.json")

## 2- run a core flooding