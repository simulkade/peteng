"""

"""
function corey_gas_oil()

end # corey_gas_oil




"""
brooks_corey_gas_oil(labda=4.17, sor=0.2, pce_go=5000)

Original Brooks and Corey (1964)

Input parameters:

    + labda: sorting factor (pore-size distribution)
    + sor: residual oil saturation
    + pce_go: gas-oil pore entry capillary pressure
"""
function brooks_corey_gas_oil(labda=4.17, sor=0.2, pce_go=5000)
  kro(so)=((so-sor)/(1-sor)).^((2+3*labda)/labda)
  krg(so)=((1-so)/(1-sor)).^2.*(1-((so-sor)/(1-sor)).^((2+labda)/labda))
  pc_go(so)=pce_go*((so-sor)/(1-sor)).^(-1/labda)
  return kro, krg, pc_go
end # brooks_corey


"""

"""
function brooks_corey_general(kro_max=1.0, krg_max=1.0, ng=2.0, no=2.0)



end # brooks_corey
