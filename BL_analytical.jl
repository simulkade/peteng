"""
returns Core-type relperms for gas and liquid as a
function of the liquid saturation and their derivatives
with respect to the liquid saturation; it also returns
the stars model for the gas relperm in presence of foam.
inputs [unit] \\(default):
   + swc [-] \\(0.1)
   + sgr [-] \\(0.05)
   + krg0 [-] \\(0.9)
   + ng [-] \\(1.8)
   + krw0 [-] \\(0.2)
   + nw [-] \\(4.2)
   + fmmob [-] \\(25000)
   + epdry [-] \\(10000)
   + fmdry [-] \\(0.29)
Outputs:
   + krw
   + krg
   + dkrwdsw
   + dkrgdsw
   + krgf
   + dkrgfdsw
Usage:
   + (krw, krg, dkrwdsw, dkrgdsw, krgf, dkrgfdsw)=corey_rel_perms()
   + (krw, krg, dkrwdsw, dkrgdsw, krgf, dkrgfdsw)=corey_rel_perms(swc,
sgr, krg0, ng, krw0, nw, fmmob, epdry, fmdry)
"""
function corey_rel_perms(;swc=0.1, sgr=0.05, krg0=0.9, ng=1.8, krw0=0.2,
nw=4.2, fmmob=25000.0, epdry=10000.0, fmdry=0.29)
  sws(sw::Real)=((sw>swc)*(sw<1-sgr)*(sw-swc)/(1-sgr-swc)+(sw>=1-sgr)*1.0)
  sws(sw::Array{Float64})=((sw.>swc).*(sw.<1-sgr).*(sw-swc)/(1-sgr-swc)+(sw.>=1-sgr).*ones(size(sw)))
  krg(sw)=((sw.>=swc).*krg0.*(1-sws(sw)).^ng+(sw.<swc).*(1+(krg0-1)/swc*sw))
  krw(sw)=((sw.<=1-sgr).*krw0.*sws(sw).^nw+(sw.>1-sgr).*(-(1-krw0)/sgr.*(1.0-sw)+1.0))
  fm(sw)=(1+fmmob*(0.5+atan(epdry.*(sw-fmdry))/π))
  krgf(sw)=(krg(sw)./fm(sw))
  dkrwdsw(sw)=(nw*krw0*(1/(1-sgr-swc))*sws(sw).^(nw-1))
  dkrgdsw(sw)=((-krg0*ng*(1-sws(sw)).^(ng-1))/(-swc-sgr+1))
  dfmdsw(sw)=(((epdry*fmmob)./(π*(epdry^2*(sw-fmdry).^2+1))))
  dkrgfdsw(sw)=((dkrgdsw(sw).*fm(sw)-dfmdsw(sw).*krg(sw))./fm(sw).^2)

  return krw, krg, dkrwdsw, dkrgdsw, krgf, dkrgfdsw
end


"""
fractional flow formulation and solution of a WAG process
for Corey-type gas and water relative permeabilities on a
Cartesian coordinate (not radial)
inputs [unit] \\(default):
   + muw [Pa.s] \\(1e-3)
   + mug [Pa.s] \\(2e-5)
   + ut [m/s] \\(1e-5)
   + phi [-] \\(0.2)
   + k [m^2] \\(1e-12)
   + swc [-] \\(0.1)
   + sgr [-] \\(0.05)
   + krg0 [-] \\(0.9)
   + ng [-] \\(1.8)
   + krw0 [-] \\(0.2)
   + nw [-] \\(4.2)
   + sw0 [-] \\(1.0)
   + sw_inj [-] \\(0.0)
   + L [m] \\(domain length = 1.0)
   + pv_inj [-] \\(injected pore volumes = 5)
Outputs:
   + xt_shock (x/t shock position)
   + sw_shock
   + xt_prf
   + sw_prf
   + t [s] \\time
   + p_inj [Pa] \\(injection pressure history)
Usage:
   using roots, Dierckx
   + (xt_shock, sw_shock, xt_prf, sw_prf, t, p_inj)=frac_flow_wag()
   + (xt_shock, sw_shock, xt_prf, sw_prf, t, p_inj)=frac_flow_wag(muw, mug, ut, phi,
k, swc, sgr, krg0, ng, krw0, nw, sw0, sw_inj, L, pv_inj)
"""
function frac_flow_wag(;muw=1e-3, mug=2e-5, ut=1e-5, phi=0.2,
  k=1e-12, swc=0.1, sgr=0.05, krg0=0.9, ng=1.8, krw0=0.2,
  nw=4.2, sw0=1.0, sw_inj=0.0, L=1, pv_inj=5)

  sws(sw::Real)=((sw>swc)*(sw<1-sgr)*(sw-swc)/(1-sgr-swc)+(sw>=1-sgr)*1.0)
  sws(sw::Array{Float64})=((sw.>swc).*(sw.<1-sgr).*(sw-swc)/(1-sgr-swc)+(sw.>=1-sgr).*ones(size(sw)))
  kr(sw)=((sw.>=swc).*krg0.*(1-sws(sw)).^ng+(sw.<swc).*(1+(krg0-1)/swc*sw))
  krw(sw)=((sw.<=1-sgr).*krw0.*sws(sw).^nw+(sw.>1-sgr).*(-(1-krw0)/sgr.*(1.0-sw)+1.0))
  krg(sw)=kr(sw)
  dkrwdsw(sw)=(nw*krw0*(1/(1-sgr-swc))*sws(sw).^(nw-1))
  dkrdsw(sw)=((-krg0*ng*(1-sws(sw)).^(ng-1))/(-swc-sgr+1))
  dkrgdsw(sw)=dkrdsw(sw)
  fw(sw)=((krw(sw)/muw)./(krw(sw)/muw+krg(sw)/mug))
  dfw(sw)=((dkrwdsw(sw)/muw.*(krw(sw)/muw+krg(sw)/mug)-
      (dkrwdsw(sw)/muw+dkrgdsw(sw)/mug).*krw(sw)/muw)./
      (krg(sw)/mug+krw(sw)/muw).^2)
    eps1=1e-3
    dfw_num(sw) = (fw(sw+eps1)-fw(sw-eps1))/(2*eps1)
    ftot(sw)=krg(sw)/mug+krw(sw)/muw

# solve the nl equation to find the shock front saturation
    f_shock(sw)=(dfw(sw)-(fw(sw0)-fw(sw))./(sw0-sw))
    sw_shock = fzero(f_shock, [swc+eps(),1-sgr-eps()])
    s=collect(linspace(0.0,1.0,100))
    s1 = collect(linspace(sw_inj, sw_shock, 500))
    xt_s1 = ut/phi*dfw(s1)
    xt_s = ut/phi*dfw(s)
    xt_shock = ut/phi*dfw(sw_shock)
    xt_prf=[xt_s1; xt_shock; xt_shock+eps(); 2*xt_shock]
    sw_prf=[s1; sw_shock; sw0; sw0]

# Go through the data first
    i=1
    while(true)
      if (i+1)==length(xt_prf)
        break
      elseif xt_prf[i]==xt_prf[i+1]
        deleteat!(xt_prf, i+1)
        deleteat!(sw_prf, i+1)
      else
        i+=1
      end
    end

# find the injection pressure history
    x = collect(linspace(0,L,1000))
    sw_int = Spline1D([xt_prf; L/eps()], [sw_prf; sw0], k=1);
    t_inj=pv_inj*phi*L/ut
    t = collect(linspace(eps(),t_inj, 100)) # [s] time
    p_inj = zeros(length(t))
    p_inj[1]=trapz(x, ut./(k*(krg(sw0*ones(size(x)))/mug+krw(sw0*ones(size(x)))/muw)))
    for i in 2:length(t)
        xt = x/t[i]
        p_inj[i] = trapz(x, ut./(k*(krg(sw_int(xt))/mug+krw(sw_int(xt))/muw)))
    end

# Return the results
    return xt_shock, sw_shock, xt_prf, sw_prf, t, p_inj
end


"""
Foam flow in porous media: a BL formulation of SAG process.
CMG Stars foam model is used for the gas mobility in
presence of foam, and Corey-type relperms are used for
the gas and water. The equations are solved on a
Cartesian coordinate in one dimension.
inputs [unit] \\(default):
   + fmmob [-] \\(25000)
   + epdry [-] \\(10000)
   + fmdry [-] \\(0.29)
   + muw [Pa.s] \\(=1e-3)
   + mug [Pa.s] \\(2e-5)
   + ut [m/s] \\(1e-5)
   + phi [-] \\(0.2)
   + k [m^2] \\(1e-12)
   + swc [-] \\(0.1)
   + sgr [-] \\(0.05)
   + krg0 [-] \\(0.9)
   + ng [-] \\(1.8)
   + krw0 [-] \\(0.2)
   + nw [-] \\(4.2)
   + sw0 [-] \\(1.0)
   + sw_inj [-] \\(0.0)
Outputs:
   + xt_shock (x/t shock position)
   + sw_shock
   + xt_prf
   + sw_prf
Usage:
   using roots
   + (xt_shock, sw_shock, xt_prf, sw_prf)=frac_flow_sag()
   + (xt_shock, sw_shock, xt_prf, sw_prf)=frac_flow_sag(fmmob, epdry, fmdry, muw, mug, ut, phi,
k, swc, sgr, krg0, ng, krw0, nw, sw0, sw_inj)
"""
function frac_flow_sag(;fmmob=25000.0, epdry=10000.0, fmdry=0.29,
  muw=1e-3, mug=2e-5, ut=1e-5, phi=0.2,
  k=1e-12, swc=0.1, sgr=0.05, krg0=0.9, ng=1.8, krw0=0.2,
  nw=4.2, sw0=1.0, sw_inj=0.0, L=1, pv_inj=5)

    sws(sw::Real)=((sw>swc)*(sw<1-sgr)*(sw-swc)/(1-sgr-swc)+(sw>=1-sgr)*1.0)
    sws(sw::Array{Float64})=((sw.>swc).*(sw.<1-sgr).*(sw-swc)/(1-sgr-swc)+(sw.>=1-sgr).*ones(size(sw)))
    kr(sw)=((sw.>=swc).*krg0.*(1-sws(sw)).^ng+(sw.<swc).*(1+(krg0-1)/swc*sw))
    krw(sw)=((sw.<=1-sgr).*krw0.*sws(sw).^nw+(sw.>1-sgr).*(-(1-krw0)/sgr.*(1.0-sw)+1.0))
    fm(sw)=(1+fmmob*(0.5+atan(epdry.*(sw-fmdry))/π))
    krg(sw)=(kr(sw)./fm(sw))
    dkrwdsw(sw)=(nw*krw0*(1/(1-sgr-swc))*sws(sw).^(nw-1))
    dkrdsw(sw)=((-krg0*ng*(1-sws(sw)).^(ng-1))/(-swc-sgr+1))
    dfmdsw(sw)=(((epdry*fmmob)./(π*(epdry^2*(sw-fmdry).^2+1))))
    dkrgdsw(sw)=((dkrdsw(sw).*fm(sw)-dfmdsw(sw).*kr(sw))./fm(sw).^2)
    fw(sw)=((krw(sw)/muw)./(krw(sw)/muw+krg(sw)/mug))
    dfw(sw)=((dkrwdsw(sw)/muw.*(krw(sw)/muw+krg(sw)/mug)-
      (dkrwdsw(sw)/muw+dkrgdsw(sw)/mug).*krw(sw)/muw)./
      (krg(sw)/mug+krw(sw)/muw).^2)
    eps1=1e-3
    dfw_num(sw) = (fw(sw+eps1)-fw(sw-eps1))/(2*eps1)
    ftot(sw)=krg(sw)/mug+krw(sw)/muw

# solve the nl equation to find the shock front saturation
    f_shock(sw)=(dfw(sw)-(fw(sw0)-fw(sw))./(sw0-sw))
    sw_shock = fzero(f_shock, [swc+eps(),1-sgr-eps()])
    s=collect(linspace(0.0,1.0,100))
    s1 = collect([linspace(sw_inj, swc, 100); linspace(swc+eps(), sw_shock, 1000)])
    xt_s1 = ut/phi*dfw(s1)
    xt_s = ut/phi*dfw(s)
    xt_shock = ut/phi*dfw(sw_shock)
    xt_prf=[xt_s1; xt_shock; xt_shock+eps(); 2*xt_shock]
    sw_prf=[s1; sw_shock; sw0; sw0]

# Go through the data first
    i=1
    while(true)
      if (i+1)==length(xt_prf)
        break
      elseif xt_prf[i]==xt_prf[i+1]
        deleteat!(xt_prf, i+1)
        deleteat!(sw_prf, i+1)
      else
        i+=1
      end
    end

# find the injection pressure history
    x = collect(linspace(0,L,1000))
    sw_int = Spline1D([xt_prf; L/eps()], [sw_prf; sw0], k=1);
    t_inj=pv_inj*phi*L/ut
    t = collect(linspace(eps(),t_inj, 100)) # [s] time
    p_inj = zeros(length(t))
    p_inj[1]=trapz(x, ut./(k*(krg(sw0*ones(size(x)))/mug+krw(sw0*ones(size(x)))/muw)))
    for i in 2:length(t)
        xt = x/t[i]
        p_inj[i] = trapz(x, ut./(k*(krg(sw_int(xt))/mug+krw(sw_int(xt))/muw)))
    end

    return xt_shock, sw_shock, xt_prf, sw_prf, t, p_inj
end

function trapz(x,y)
# copied from https://github.com/hwborchers/NumericalMath.jl/blob/master/src/integrate.jl
  n=length(x)
  if (length(y) != n)
    error("Vectors 'x', 'y' must be of same length")
  end
  r=0
  for i in 2:n
    r+=(x[i] - x[i-1])*(y[i] + y[i-1])
  end
  return r/2.0
end
