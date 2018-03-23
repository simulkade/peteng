# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 13:11:42 2018

@author: ehsan
"""

#import relative_permeability as relperm
#import matplotlib.pyplot as plt
#import numpy as np
#
#swc=0.1
#sor=0.05
#kro0=0.9
#no=2.0
#krw0=0.4
#nw=2.0
#  
#kro_new = lambda sw: relperm.kro(sw, kro0, sor, swc, no)
#krw_new = lambda sw: relperm.krw(sw, krw0, sor, swc, nw)
#dkrwdsw_new = lambda sw: relperm.dkrwdsw(sw, krw0, sor, swc, nw)
#dkrodsw_new = lambda sw: relperm.dkrodsw(sw, kro0, sor, swc, no)
#sw = np.linspace(0.0, 1.0, 100)
#plt.plot(sw, krw_new(sw), sw, kro_new(sw))
#plt.show()

import fractional_flow as ff
xt_shock, sw_shock, xt_prf, sw_prf, t, p_inj, R_oil = ff.frac_flow_wf()

import matplotlib.pyplot as plt

plt.figure()
plt.plot(t, R_oil)
plt.show()

plt.figure()
plt.plot(xt_prf, sw_prf)
plt.show()

plt.figure()
plt.plot(t, p_inj)
plt.show()