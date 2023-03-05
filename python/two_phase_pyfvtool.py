# def two_phase_flow(muw=1e-3, muo=2e-3, ut=1e-5, phi=0.2,
#                    k=1e-12, swc=0.1, sor=0.05, kro0=0.9, no=2.0, krw0=0.4,
#                    nw=2.0, sw0=0.0, sw_inj=1.0, L=1.0, pv_inj=5.0):
# Coupled nonlinear PDE's
# Buckley Leverett equation
# dependent variables: pressure and water saturation
# Prepared for educational purposes by ** AAE **
# works fine, timestepping can be improved (ODE solver?)
# Written by Ali A. Eftekhari
# Last checked: June 2021

from pyfvtool import *
import matplotlib.pyplot as plt
import numpy as np

class Fluids:
    def __init__(self, mu_water=0.001, mu_oil=0.003):
        self.water_viscosity = mu_water
        self.oil_viscosity = mu_oil


class RelativePermeability:
    def __init__(self, swc=0.1, sor=0.05, kro0=0.9, no=2.0,
                 krw0=0.4, nw=2.0):
        self.kro0 = kro0
        self.krw0 = krw0
        self.no = no
        self.nw = nw
        self.swc = swc
        self.sor = sor

    def kro(self, sw):
        kro0 = self.kro0
        sor = self.sor
        swc = self.swc
        no = self.no
        res = ((swc <= sw) & (sw <= 1-sor))*kro0*((1-sw-sor)/(1-sor-swc))**no \
        + ((0.0 < sw) & (sw < swc))*(1+(kro0-1)/swc*sw) \
        + (sw > 1-sor)*0.0 \
            + (sw <= 0.0)*1.0
        return res

    def krw(self, sw):
        krw0, sor, swc, nw = self.krw0, self.sor, self.swc, self.nw
        res = ((swc <= sw) & (sw <= 1-sor))*krw0*((sw-swc)/(1-sor-swc))**nw \
        + ((1-sor < sw) & (sw < 1.0))*(-(1-krw0)/sor*(1.0-sw)+1.0) \
        + (sw <= swc)*0.0 \
            + (sw >= 1.0)*1.0
        return res

    def dkrodsw(self, sw):
        krw0, sor, swc, nw = self.krw0, self.sor, self.swc, self.nw
        res = ((swc <= sw) & (sw <= 1-sor))*(-no*kro0/(1-sor-swc)*((1-sw-sor)/(1-sor-swc))**(no-1)) \
        + ((0.0 < sw) & (sw < swc))*(kro0-1)/swc \
            + ((sw > 1-sor) | (sw <= 0.0))*0.0
        return res

    def dkrwdsw(self, sw):
        krw0, sor, swc, nw = self.krw0, self.sor, self.swc, self.nw
        res = ((swc <= sw) & (sw <= 1-sor))*nw*krw0/(1-sor-swc)*((sw-swc)/(1-sor-swc))**(nw-1) \
        + ((1-sor < sw) & (sw < 1.0))*(1-krw0)/sor \
            + ((sw < swc) | (sw >= 1.0))*0.0
        return res
    
    def visualise(self):
       plt.figure()
       sw_ = np.linspace(0.0, 1.0, 50)
       plt.plot(sw_, self.krw(sw_), label="Water")
       plt.plot(sw_, self.kro(sw_), label="Oil")
       plt.xlabel("Water saturation")
       plt.ylabel("Relative permeability")
       plt.legend()


class Resevoir:
    def __init__(self, rel_perm: RelativePermeability, fluids: Fluids, 
                 porosity=0.2, permeability=0.01e-12,
                 sw_init = 0.2, pressure_init = 100e5):
        self.porosity = porosity
        self.permeability = permeability
        self.rel_perm = rel_perm
        self.fluids = fluids
        self.initial_sw = np.maximum(rel_perm.swc, sw_init)
        self.initial_p = pressure_init

class ReservoirModel1D:
    def __init__(self, Nx: int, length=300.0,
                 injection_velocity = 1e-5,
                 dp_allowed = 100, dsw_allowed = 0.05,
                 eps_p = 1e-5, eps_sw = 1e-5):
        m = createMesh1D(Nx, length)
        self.domain = m
        BCp = createBC(m)  # Neumann BC for pressure
        BCs = createBC(m)  # Neumann BC for saturation


        
        self.pressure_bc = BCp
        self.saturation_bc = BCs
        # self.BC = ...
        # self.IC = ...
        # self.


# define the geometry
Nx = 100  # number of cells in x direction
W = 300  # [m] length of the domain in x direction
m = createMesh1D(Nx, W)  # creates a 1D mesh
# define the physical parametrs
krw0 = 1.0
kro0 = 0.76
nw = 2.4
no = 2.0
sor = 0.12
swc = 0.09
sws = @(sw)((sw > swc).*(sw < 1-sor).*(sw-swc)/(1-sor-swc)+(sw >= 1-sor).*ones(size(sw)))
kro = @(sw)((sw >= swc).*kro0.*(1-sws(sw)). ^ no+(sw < swc).*(1+(kro0-1)/swc*sw))
krw = @(sw)((sw <= 1-sor).*krw0.*sws(sw). ^ nw+(sw > 1-sor).*(-(1-krw0)/sor.*(1.0-sw)+1.0))
dkrwdsw = @(sw)((sw <= 1-sor).*nw.*krw0.*(1/(1-sor-swc)).*sws(sw). ^ (nw-1)+(sw > 1-sor)*((1-krw0)/sor))
dkrodsw = @(sw)((sw >= swc).*(-kro0*no*(1-sws(sw)). ^ (no-1))/(-swc-sor+1)+(sw < swc).*((kro0-1)/swc))
p0 = 100e5  # [bar] pressure
pin = 150e5  # [bar] injection pressure at the left boundary
u_in = 1.0/(24*3600)  # [m/s] equal to 1 m/day
sw0 = swc+0.1  # initial water saturation
sw_in = 1
mu_oil = 2e-3  # [Pa.s] oil viscosity
mu_water = 1e-3  # [Pa.s] water viscosity
# reservoir
k0 = 2e-12  # [m^2] average reservoir permeability
phi0 = 0.2  # average porosity
clx = 1.2
cly = 0.2
V_dp = 0.7  # Dykstra-Parsons coef.
perm_val = k0  # field2d(Nx,Ny,k0,V_dp,clx,cly)
k = createCellVariable(m, perm_val)
phi = createCellVariable(m, phi0)
lw = geometricMean(k)/mu_water
lo = geometricMean(k)/mu_oil
Lw = @(sw)(krw(sw))
Lo = @(sw)(k/mu_oil*kro(sw))
dLwdsw = @(sw)(k/mu_water*dkrwdsw(sw))
dLodsw = @(sw)(k/mu_oil*dkrodsw(sw))
# Define the boundaries
BCp = createBC(m)  # Neumann BC for pressure
BCs = createBC(m)  # Neumann BC for saturation
# left boundary pressure gradient
BCp.left.a[:] = (krw(sw_in)*lw.xvalue(1, :)+kro(sw_in)*lo.xvalue(1, : )) BCp.left.b(: ) = 0 BCp.left.c(: ) = -u_in
# change the right boandary to constant pressure (Dirichlet)
# BCp.left.a[:]=0 BCp.left.b[:]=1 BCp.left.c[:]=pin
BCp.right.a[:] = 0 BCp.right.b(: ) = 1 BCp.right.c(: ) = p0
# change the left boundary to constant saturation (Dirichlet)
BCs.left.a[:] = 0 BCs.left.b(: ) = 1 BCs.left.c[:] = 1
# define the time step and solver properties
# dt = 1000 # [s] time step
dt = (W/Nx)/u_in/10  # [s]
t_end = 1000*dt  # [s] final time
eps_p = 1e-5  # pressure accuracy
eps_sw = 1e-5  # saturation accuracy
# define the variables
sw_old = createCellVariable(m, sw0, BCs)
p_old = createCellVariable(m, p0, BCp)
sw = sw_old
p = p_old
uw = -gradientTerm(p_old)  # an estimation of the water velocity
# start the main loop
# generate intial pressure profile (necessary to initialize the fully
# implicit solver)
dp_alwd = 100.0  # Pa
dsw_alwd = 0.05
t = 0
# fprintf(1, 'progress (##):  ')
while (t < t_end)
 error_p = 1e5
  error_sw = 1e5
   # Implicit loop
   loop_count = 0
    while ((error_p > eps_p) | | (error_sw > eps_sw))
      loop_count = loop_count+1
       if loop_count > 10
         break
        end
        # calculate parameters
        pgrad = gradientTerm(p)
        sw_face = upwindMean(sw, -pgrad)  # average value of water saturation
        labdao = lo.*funceval(kro, sw_face)
        labdaw = lw.*funceval(krw, sw_face)
        dlabdaodsw = lo.*funceval(dkrodsw, sw_face)
        dlabdawdsw = lw.*funceval(dkrwdsw, sw_face)
        labda = labdao+labdaw
        dlabdadsw = dlabdaodsw+dlabdawdsw
        # compute [Jacobian] matrices
        Mdiffp1 = diffusionTerm(-labda)
        Mdiffp2 = diffusionTerm(-labdaw)
        Mconvsw1 = convectionUpwindTerm(-dlabdadsw.*pgrad)
        Mconvsw2 = convectionUpwindTerm(-dlabdawdsw.*pgrad)
        [Mtranssw2, RHStrans2] = transientTerm(sw_old, dt, phi)
        # Compute RHS values
        RHS1 = divergenceTerm(-dlabdadsw.*sw_face.*pgrad)
        RHS2 = divergenceTerm(-dlabdawdsw.*sw_face.*pgrad)
        # include boundary conditions
        [Mbcp, RHSbcp] = boundaryCondition(BCp)
        [Mbcsw, RHSbcsw] = boundaryCondition(BCs)
        # Couple the equations BC goes into the block on the main diagonal
        M = [Mdiffp1+Mbcp Mconvsw1 Mdiffp2 Mconvsw2+Mtranssw2+Mbcsw]
        RHS = [RHS1+RHSbcp RHS2+RHStrans2+RHSbcsw]
        # solve the linear system of equations
        x = M\RHS
        # x = agmg(M, RHS, [], 1e-10, 500, [], [p.value[:] sw.value[:]])
        # separate the variables from the solution
        p_new = reshapeCell(m, full(x(1: (Nx+2))))
        sw_new = reshapeCell(m, full(x((Nx+2)+1: end)))
        # calculate error values
        error_p = max(abs((p_new(: )-p.value(: ))./p_new(: )))
        error_sw = max(abs(sw_new[:]-sw.value(: )))
        # assign new values of p and sw
        p.value = p_new
        sw.value = sw_new
    end
    if loop_count > 10
      p = p_old
       sw = sw_old
        dt = dt/5
        continue
    end
    dsw = max(abs(sw_new(: )-sw_old.value[:])./sw_new[:])
    t = t+dt
    fprintf(1, '\b\b\b\b\b\b\b\b\b\b\b\b\b\bProgress: #d ##', floor(t/t_end*100))
    dt = min([dt*(dsw_alwd/dsw), 2*dt, t_end-t])
    p_old = p
    sw_old = sw
    figure(1)visualizeCells(sw) drawnow
end
fprintf('\n')
