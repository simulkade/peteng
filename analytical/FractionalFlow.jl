module FractionalFlow

using Roots, Dierckx

# types:
struct CoreyRelativePermeability
    krw::Real
    kro::Real
    swc::Real
    sor::Real
    nw::Real
    no::Real
end

"""
"""
function oil_water_rel_perm(;krw=0.4, kro=0.9, 
    swc=0.15, sor=0.2, nw=2.0, no = 2.0)
    CoreyRelativePermeability(krw, kro, swc, sor, nw, no)
end

end # module