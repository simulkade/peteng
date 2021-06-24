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
    res=-no*kro0/(1-sor-swc)*((1-sw-sor)/(1-sor-swc))^(no-1)
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
    res=nw*krw0/(1-sor-swc)*((sw-swc)/(1-sor-swc))^(nw-1)
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
  pc0=1.0e7
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
  pc0=1.0e7
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

# all rel-perm functions for array Inputs -------------------------------------
function sws(sw::Array{Float64}, sor::Array{Float64}, swc::Array{Float64})
  res=zeros(size(sw))
  for i in eachindex(sw)
    res[i]=sws(sw[i], sor[i], swc[i])
  end
  return res
end

function kro(sw::Array{Float64}, kro0::Array{Float64}, sor::Array{Float64}, swc::Array{Float64}, no::Array{Float64})
  res=zeros(size(sw))
  for i in eachindex(sw)
    res[i]=kro(sw[i], kro0[i], sor[i], swc[i], no[i])
  end
  return res
end

function dkrodsw(sw::Array{Float64}, kro0::Array{Float64}, sor::Array{Float64}, swc::Array{Float64}, no::Array{Float64})
  res=zeros(size(sw))
  for i in eachindex(sw)
    res[i]=dkrodsw(sw[i], kro0[i], sor[i], swc[i], no[i])
  end
  return res
end

function krw(sw::Array{Float64}, krw0::Array{Float64}, sor::Array{Float64}, swc::Array{Float64}, nw::Array{Float64})
  res=zeros(size(sw))
  for i in eachindex(sw)
    res[i]=krw(sw[i], krw0[i], sor[i], swc[i], nw[i])
  end
  return res
end

function dkrwdsw(sw::Array{Float64}, krw0::Array{Float64}, sor::Array{Float64}, swc::Array{Float64}, nw::Array{Float64})
  res=zeros(size(sw))
  for i in eachindex(sw)
    res[i]=dkrwdsw(sw[i], krw0[i], sor[i], swc[i], nw[i])
  end
  return res
end

function pc_drain(sw::Array{Float64}, pce::Array{Float64}, swc::Array{Float64}; labda::Real=2.4)
  res=zeros(size(sw))
  for i in eachindex(sw)
    res[i]=pc_drain(sw[i], pce[i], swc[i], labda=labda)
  end
  return res
end

function dpc_drain(sw::Array{Float64}, pce::Array{Float64}, swc::Array{Float64}; labda::Real=2.4)
  res=zeros(size(sw))
  for i in eachindex(sw)
    res[i]=dpc_drain(sw[i], pce[i], swc[i], labda=labda)
  end
  return res
end

function pc_imb(sw::Array{Float64}, pce::Array{Float64}, swc::Array{Float64}, sor::Array{Float64},
  teta::Array{Float64}; labda::Real=2.4, b::Real=0.6)
  res=zeros(size(sw))
  for i in eachindex(sw)
    res[i]=pc_imb(sw[i], pce[i], swc[i], sor[i], teta=teta[i], labda=labda, b=b)
  end
  return res
end

function dpc_imb(sw::Array{Float64}, pce::Array{Float64}, swc::Array{Float64}, sor::Array{Float64},
  teta::Array{Float64}; labda::Real=2.4, b::Real=0.6)
  res=zeros(size(sw))
  for i in eachindex(sw)
    res[i]=dpc_imb(sw[i], pce[i], swc[i], sor[i], teta=teta[i], labda=labda, b=b)
  end
  return res
end

# sws=@(sw, sor, swc)((sw>swc).*(sw<1-sor).*(sw-swc)./(1-sor-swc)+(sw>=1-sor))
# kro=@(sw, kro0, sor, swc)((sw>=swc).*kro0.*(1-sws(sw, sor, swc)).^no+(sw<swc).*(1+(kro0-1)./swc.*sw))
# krw=@(sw, krw0, sor, swc)((sw<=1-sor).*krw0.*sws(sw, sor, swc).^nw+(sw>1-sor).*(-(1-krw0)./sor.*(1.0-sw)+1.0))
# dkrwdsw=@(sw, krw0, sor, swc)((sw<=1-sor).*nw.*krw0.*(1./(1-sor-swc)).*sws(sw, sor, swc).^(nw-1)+(sw>1-sor).*((1-krw0)./sor))
# dkrodsw=@(sw, kro0, sor, swc)((sw>=swc).*(-kro0.*no.*(1-sws(sw, sor, swc)).^(no-1))./(-swc-sor+1)+(sw<swc).*((kro0-1)./swc))
