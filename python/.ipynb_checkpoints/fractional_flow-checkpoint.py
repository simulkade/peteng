import relative_permeability as relperm
import numpy as np
import scipy.optimize as opt
from scipy.interpolate import interp1d

def frac_flow_wf(muw=1e-3, muo=2e-3, ut=1e-5, phi=0.2, \
  k=1e-12, swc=0.1, sor=0.05, kro0=0.9, no=2.0, krw0=0.4, \
  nw=2.0, sw0=0.0, sw_inj=1.0, L=1.0, pv_inj=5.0):
  
    # sws(sw::Real)=((sw>swc)*(sw<1-sor)*(sw-swc)/(1-sor-swc)+(sw>=1-sor)*1.0)
    # sws(sw::Array{Float64})=((sw.>swc).*(sw.<1-sor).*(sw-swc)/(1-sor-swc)+(sw.>=1-sor).*ones(size(sw)))
    # kro_new(sw)=((sw.>=swc).*kro0.*(1-sws(sw)).^no+(sw.<swc).*(1+(kro0-1)/swc*sw))
    # krw_new(sw)=((sw.<=1-sor).*krw0.*sws(sw).^nw+(sw.>1-sor).*(-(1-krw0)/sor.*(1.0-sw)+1.0))
    # dkrwdsw_new(sw)=(nw*krw0*(1/(1-sor-swc))*sws(sw).^(nw-1))
    # dkrodsw_new(sw)=((-kro0*no*(1-sws(sw)).^(no-1))/(-swc-sor+1))
    kro_new = lambda sw: relperm.kro(sw, kro0, sor, swc, no)
    krw_new = lambda sw: relperm.krw(sw, krw0, sor, swc, nw)
    dkrwdsw_new = lambda sw: relperm.dkrwdsw(sw, krw0, sor, swc, nw)
    dkrodsw_new = lambda sw: relperm.dkrodsw(sw, kro0, sor, swc, no)
    fw = lambda sw: ((krw_new(sw)/muw)/(krw_new(sw)/muw+kro_new(sw)/muo))
    dfw = lambda sw: ((dkrwdsw_new(sw)/muw*(krw_new(sw)/muw+kro_new(sw)/muo)- \
    (dkrwdsw_new(sw)/muw+dkrodsw_new(sw)/muo)*krw_new(sw)/muw)/ \
    (kro_new(sw)/muo+krw_new(sw)/muw)**2)
    # eps1=1e-3
    # dfw_num(sw) = (fw(sw+eps1)-fw(sw-eps1))/(2*eps1)
    # ftot(sw)=kro_new(sw)/muo+krw_new(sw)/muw

    # solve the nl equation to find the shock front saturation
    eps = 1e-15
    f_shock = lambda sw: (dfw(sw)-(fw(sw)-fw(sw0))/(sw-sw0))
    sw_tmp = np.linspace(swc+eps, 1-sor-eps, 1000)
    sw_shock = sw_tmp[np.abs(f_shock(sw_tmp)).argmin()]
    try:
        if f_shock(swc+eps)*f_shock(1-sor-eps)>0:
            sw_shock = opt.newton(f_shock, sw_shock) #  ftol = 1e-10, xtol = 1e-7)
        else:
            sw_shock = opt.brentq(f_shock, swc+eps, 1-sor-eps)
    except:
        print('shock front saturation is estimated: $sw_shock, error is $(f_shock(sw_shock))')
    s = np.linspace(0.0, 1.0, 100)
    s1 = np.linspace(sw_inj, sw_shock-eps, 1000)
    xt_s1 = ut/phi*dfw(s1)
    xt_shock = ut/phi*dfw(sw_shock)
    xt_prf=np.concatenate((xt_s1, [xt_shock, xt_shock+eps, 2*xt_shock]))
    sw_prf=np.concatenate((s1, [sw_shock, sw0, sw0]))
    # println(xt_prf)
    # print(sw_shock)
    # Go through the data first
    i=1
    xt_prf, indices = np.unique(xt_prf, return_index = True)
    sw_prf = sw_prf[indices]
    
    # xt_prf_slim = unique(xt_prf)
    # sw_prf = sw_prf[indexin(xt_prf_slim, xt_prf)]
    # xt_prf = xt_prf_slim

    # find the injection pressure history
    x = np.linspace(0,L,1000)
    # sw_int = Spline1D([xt_prf; L/eps()], [sw_prf; sw0], k=1)
    # println(xt_prf)
    # println(sw_prf)
    sw_int = interp1d(xt_prf, sw_prf, kind='linear', fill_value='extrapolate')
    t_inj=pv_inj*phi*L/ut
    t = np.linspace(0.0,t_inj, 200) # [s] time
    p_inj = np.zeros(t.size)
    R_oil= np.zeros(t.size)
    p_inj[0]=np.trapz(ut/(k*(kro_new(sw0*np.ones(np.size(x)))/muo+krw_new(sw0*np.ones(np.size(x)))/muw)), x=x)
    for i in range(1,t.size):
        xt = x/t[i]
        p_inj[i] = np.trapz(ut/(k*(kro_new(sw_int(xt))/muo+krw_new(sw_int(xt))/muw)), x=x)
        R_oil[i]=1.0-np.trapz(1.0-sw_int(xt), x = x/L)/(1-sw0)

    # Return the results
    return (xt_shock, sw_shock, xt_prf, sw_prf, t, p_inj, R_oil)
