"""
normalized water saturation
"""
function sws(sw::Real, sor::Real, swc::Real)
  if sw>swc && sw<1-sor
    res=(sw-swc)/(1-sor-swc)
  elseif sw>=1-sor
    res=1.0
  elseif sw<=swc
    res=0.0
  end
  return res
end

"""
Oil relative permeability functions for oil-water system
"""
function kro(sw::Real, kro0::Real, sor::Real, swc::Real, no::Real)
  if swc<=sw<=1-sor
    res=kro0*((1-sw-sor)/(1-sor-swc))^no
  elseif 0.0<sw<swc
    res=1+(kro0-1)/swc*sw
  elseif sw>1-sor
    res=0.0
  elseif sw<=0.0
    res=1.0
  end
  return res
end

function dkrodsw(sw::Real, kro0::Real, sor::Real, swc::Real, no::Real)
  if swc<=sw<=1-sor
    res=-kro0/(1-sor-swc)*((1-sw-sor)/(1-sor-swc))^(no-1)
  elseif 0.0<sw<swc
    res=(kro0-1)/swc
  elseif sw>1-sor || sw<=0.0
    res=0.0
  end
  return res
end

"""
water rel-perm for oil-water system
"""
function krw(sw::Real, krw0::Real, sor::Real, swc::Real, nw::Real)
  if swc<=sw<=1-sor
    res=krw0*((sw-swc)/(1-sor-swc))^nw
  elseif 1-sor<sw<1.0
    res=-(1-krw0)/sor*(1.0-sw)+1.0
  elseif sw<swc
    res=0.0
  elseif sw>=1.0
    res=1.0
  end
  return res
end

function dkrwdsw(sw::Real, krw0::Real, sor::Real, swc::Real, nw::Real)
  if swc<=sw<=1-sor
    res=krw0/(1-sor-swc)*((sw-swc)/(1-sor-swc))^(nw-1)
  elseif 1-sor<sw<1.0
    res=(1-krw0)/sor
  elseif sw<swc || sw>=1.0
    res=0.0
  end
  return res
end

"""
capillary pressure curve for drainage
"""
function pc_drain(sw::Real, pce::Real, swc::Real; labda::Real=2.4)
  pc0=1.0e9
  sw0=swc+(1-labda*log(pc0/pce)+sqrt((-1+labda*log(pc0/pce))^2+4*swc/(1-swc)))/2*(1-swc)
  if sw>sw0
    res=pce*((sw-swc)/(1-swc))^(-1.0/labda)
  elseif 0.0<=sw<=sw0
    pcs=pce*((sw0-swc)/(1-swc))^(-1.0/labda)
    res=exp((log(pcs)-log(pc0))/sw0*(sw-sw0)+log(pcs))
  elseif sw<0
    res=pc0
  else
    res=0.0
  end
end

function dpc_drain(sw::Real, pce::Real, swc::Real; labda::Real=2.4)
  pc0=1.0e9
  sw0=swc+(1-labda*log(pc0/pce)+sqrt((-1+labda*log(pc0/pce))^2+4*swc/(1-swc)))/2*(1-swc)
  if sw>sw0
    res=-1.0/((1-swc)*labda)*pce*((sw-swc)/(1-swc))^(-1.0/labda-1)
  elseif 0.0<=sw<=sw0
    res=-1.0/((1-swc)*labda)*pce*((sw0-swc)/(1-swc))^(-1.0/labda-1)
  else
    res=0.0
  end
end

function pc_imb(sw::Real, pce::Real, swc::Real, sor::Real;
  teta::Real=0.785, labda::Real=2.4, b::Real=0.6)
  pc1=pc_drain(sw, pce, swc, labda=labda)
  pc2=pc_drain(1-sw, pce, sor, labda=labda)
  return (0.5*(1+cos(teta)))^b*pc1-(0.5*(1-cos(teta)))^b*pc2
end

function dpc_imb(sw::Real, pce::Real, swc::Real, sor::Real; teta::Real=0.785, labda::Real=2.4, b::Real=0.6)
  dpc1=dpc_drain(sw, pce, swc, labda=labda)
  dpc2=dpc_drain(1-sw, pce, sor, labda=labda)
  return (0.5*(1+cos(teta)))^b*dpc1+(0.5*(1-cos(teta)))^b*dpc2
end
