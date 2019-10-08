# relative permeability module
import numpy as np

# def kro(sw, kro0, sor, swc, no):
#     if swc<=sw<=1-sor:
#         res=kro0*((1-sw-sor)/(1-sor-swc))**no
#     elif 0.0<sw<swc:
#         res=1+(kro0-1)/swc*sw
#     elif sw>1-sor:
#         res=0.0
#     elif sw<=0.0:
#         res=1.0
#     return res


def kro(sw, kro0, sor, swc, no):
    res = ((swc<=sw) & (sw<=1-sor))*kro0*((1-sw-sor)/(1-sor-swc))**no \
    +((0.0<sw) & (sw<swc))*(1+(kro0-1)/swc*sw) \
    +(sw>1-sor)*0.0 \
    +(sw<=0.0)*1.0
    return res

def krw(sw, krw0, sor, swc, nw):
    res = ((swc<=sw) & (sw<=1-sor))*krw0*((sw-swc)/(1-sor-swc))**nw \
    +((1-sor<sw) & (sw<1.0))*(-(1-krw0)/sor*(1.0-sw)+1.0) \
    +(sw<=swc)*0.0 \
    +(sw>=1.0)*1.0
    return res

def dkrodsw(sw, kro0, sor, swc, no):
    res = ((swc<=sw) & (sw<=1-sor))*(-no*kro0/(1-sor-swc)*((1-sw-sor)/(1-sor-swc))**(no-1)) \
    +((0.0<sw) & (sw<swc))*(kro0-1)/swc \
    +((sw>1-sor) | (sw<=0.0))*0.0
    return res

def dkrwdsw(sw, krw0, sor, swc, nw):
    res = ((swc<=sw) & (sw<=1-sor))*nw*krw0/(1-sor-swc)*((sw-swc)/(1-sor-swc))**(nw-1) \
    +((1-sor<sw) & (sw<1.0))*(1-krw0)/sor \
    +((sw<swc) | (sw>=1.0))*0.0
    return res