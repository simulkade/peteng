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
from scipy.sparse import vstack, hstack
from scipy.sparse.linalg import spsolve

class Fluids:
    """
    This class represents a model for calculating fluid properties.

    Attributes:
        water_viscosity (float): Viscosity of water. Default is 0.001 Pa.s.
        oil_viscosity (float): Viscosity of oil. Default is 0.003 Pa.s.
    """
    def __init__(self, mu_water=0.001, mu_oil=0.003):
        self.water_viscosity = mu_water
        self.oil_viscosity = mu_oil


class RelativePermeability:
    """
    This class represents a model for calculating relative permeability values
    based on water saturation (sw) using the Corey-Brooks model.

    Attributes:
        swc (float): Connate water saturation. Default is 0.1.
        sor (float): Residual oil saturation. Default is 0.05.
        kro0 (float): Oil relative permeability at 100% water saturation. Default is 0.9.
        no (float): Corey exponent for oil relative permeability. Default is 2.0.
        krw0 (float): Water relative permeability at 100% water saturation. Default is 0.4.
        nw (float): Corey exponent for water relative permeability. Default is 2.0.
    """

    def __init__(self, swc=0.1, sor=0.05, kro0=0.9, no=2.0, krw0=0.4, nw=2.0):
        """
        Initializes a new instance of the RelativePermeability class.

        Args:
            swc (float, optional): Connate water saturation. Default is 0.1.
            sor (float, optional): Residual oil saturation. Default is 0.05.
            kro0 (float, optional): Oil relative permeability at 100% water saturation. Default is 0.9.
            no (float, optional): Corey exponent for oil relative permeability. Default is 2.0.
            krw0 (float, optional): Water relative permeability at 100% water saturation. Default is 0.4.
            nw (float, optional): Corey exponent for water relative permeability. Default is 2.0.
        """
        self.kro0 = kro0
        self.krw0 = krw0
        self.no = no
        self.nw = nw
        self.swc = swc
        self.sor = sor

    def kro(self, sw:np.ndarray) -> np.ndarray:
        """
        Calculates the oil relative permeability based on the water saturation.

        Args:
            sw (float): Water saturation.

        Returns:
            float: Oil relative permeability.
        """
        kro0 = self.kro0
        sor = self.sor
        swc = self.swc
        no = self.no
        res = np.zeros_like(sw)
        cond1 = np.logical_and(swc <= sw, sw <= 1 - sor)
        res[cond1] = kro0 * ((1 - sw[cond1] - sor) / (1 - sor - swc)) ** no
        cond2 = np.logical_and(0.0 < sw, sw < swc)
        res[cond2] = 1 + (kro0 - 1) / swc * sw[cond2]
        res[sw > 1 - sor] = 0.0
        res[sw <= 0.0] = 1.0
        return res

    def krw(self, sw:np.ndarray) -> np.ndarray:
        """
        Calculates the water relative permeability based on the water saturation.

        Args:
            sw (float): Water saturation.

        Returns:
            float: Water relative permeability.
        """
        krw0, sor, swc, nw = self.krw0, self.sor, self.swc, self.nw
        res = np.zeros_like(sw)
        cond1 = np.logical_and(swc <= sw, sw <= 1 - sor)
        res[cond1] = krw0 * ((sw[cond1] - swc) / (1 - sor - swc)) ** nw
        cond2 = np.logical_and(1 - sor < sw, sw < 1.0)
        res[cond2] = (-(1 - krw0) / sor * (1.0 - sw[cond2]) + 1.0)
        res[sw <= swc] = 0.0
        res[sw >= 1.0] = 1.0
        return res

    def dkrodsw(self, sw):
        """
        Calculates the derivative of oil relative permeability with respect to water saturation.

        Args:
            sw (float): Water saturation.

        Returns:
            float: Derivative of oil relative permeability with respect to water saturation.
        """
        kro0, sor, swc, no = self.kro0, self.sor, self.swc, self.no
        res = np.zeros_like(sw)
        cond1 = np.logical_and(swc <= sw, sw <= 1 - sor)
        res[cond1] = -no * kro0 / (1 - sor - swc) * ((1 - sw[cond1] - sor) / (1 - sor - swc)) ** (no - 1)
        cond2 = np.logical_and(0.0 < sw, sw < swc)
        res[cond2] = (kro0 - 1) / swc
        res[sw > 1 - sor] = 0.0 
        res[sw <= 0.0] = 0.0
        return res

    def dkrwdsw(self, sw):
        """
        Calculates the derivative of water relative permeability with respect to water saturation.

        Args:
            sw (float): Water saturation.

        Returns:
            float: Derivative of water relative permeability with respect to water saturation.
        """
        krw0, sor, swc, nw = self.krw0, self.sor, self.swc, self.nw
        res = np.zeros_like(sw)
        cond1 = np.logical_and(swc <= sw, sw <= 1 - sor)
        res[cond1] = nw * krw0 / (1 - sor - swc) * ((sw[cond1] - swc) / (1 - sor - swc)) ** (nw - 1)
        cond2 = np.logical_and(1 - sor < sw, sw < 1.0)
        res[cond2] = (1 - krw0) / sor
        res[sw < swc] = 0.0
        res[sw >= 1.0] = 0.0
        return res

    def visualize(self):
        """
        Visualizes the relative permeability curves for water and oil.

        Requires matplotlib.pyplot to be imported as plt.

        Returns:
            None
        """
        import matplotlib.pyplot as plt
        import numpy as np

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

class OperationalConditions:
    """
    Operational conditions for the reservoir. This includes the injection rate and 
    injection/production pressures.

    Args:
        injection_velocity (float, optional): Injection velocity in m/s. Default is 1e-5.
        injection_pressure (float, optional): Injection pressure in Pa. Default is 100e5.
        production_pressure (float, optional): Production pressure in Pa. Default is 50e5.
    """
    def __init__(self, injection_velocity: float=1e-5,               
                 injection_pressure: float=100e5,                 
                 production_pressure : float=50e5,
                 injection_sw = 1.0,
                 active_rate= True):
        self.injection_velocity = injection_velocity
        self.injection_pressure = injection_pressure
        self.production_pressure = production_pressure
        self.active_rate = active_rate
        self.injection_sw = injection_sw

class CoreModel1D:
    def __init__(self, 
                 reservoir: Resevoir,
                 operational_conditions: OperationalConditions, 
                 Nx: int = 50, 
                 length: float=6.0e-2,
                 diameter: float=2.5e-2,
                 dp_allowed = 100, dsw_allowed = 0.05,
                 eps_p = 1e-5, eps_sw = 1e-5):
        m = createMesh1D(Nx, length)
        self.pore_volume = reservoir.porosity * length * np.pi * (diameter/2)**2
        self.injection_velocity = operational_conditions.injection_velocity
        self.dp_allowed = dp_allowed
        self.dsw_allowed = dsw_allowed
        self.eps_p = eps_p
        self.eps_sw = eps_sw
        self.domain = m
        self.perm_field = createCellVariable(m, reservoir.permeability)
        self.poros_field = createCellVariable(m, reservoir.porosity)
        mu_water = createCellVariable(m, reservoir.fluids.water_viscosity)
        mu_oil = createCellVariable(m, reservoir.fluids.oil_viscosity)
        self.water_mobility_max = geometricMean(self.perm_field/mu_water)
        self.oil_mobility_max = geometricMean(self.perm_field/mu_oil)
        self.rel_perm = reservoir.rel_perm
        
        BCp = createBC(m)  # Neumann BC for pressure
        # back pressure boundary condition
        BCp.right.a[:] = 0.0
        BCp.right.b[:] = 1.0
        BCp.right.c[:] = operational_conditions.production_pressure
        # injection boundary
        BCp.left.a[:] = self.water_mobility_max[0]
        BCp.left.b[:] = 0.0
        BCp.left.c[:] = -operational_conditions.injection_velocity
        # saturation left boundary
        BCs = createBC(m)  # Neumann BC for saturation
        BCs.left.a[:] = 0.0
        BCs.left.b[:] = 1.0
        BCs.left.c[:] = operational_conditions.injection_sw

        self.pressure_bc = BCp
        self.saturation_bc = BCs
        self.initial_pressure = createCellVariable(m, reservoir.initial_p, BCp)
        self.initial_sw = createCellVariable(m, reservoir.initial_sw, BCs)
        self.pressure = self.initial_pressure
        self.saturation = self.initial_sw

    def calculate_water_mobility(self, sw_face):
        return self.water_mobility_max * faceeval(self.rel_perm.krw, sw_face)
    
    def calculate_oil_mobility(self, sw_face):
        return self.oil_mobility_max * faceeval(self.rel_perm.kro, sw_face)
    
    def calculate_water_mobility_derivative(self, sw_face):
        return self.water_mobility_max * faceeval(self.rel_perm.dkrwdsw, sw_face)
    
    def calculate_oil_mobility_derivative(self, sw_face):
        return self.oil_mobility_max * faceeval(self.rel_perm.dkrodsw, sw_face)
    
    def calculate_water_velocity(self, sw_face):
        return self.calculate_water_mobility(sw_face) * gradientTerm(self.pressure)
    
    def simulate(self, final_time = 100000.0, dt = 1000.0):
        t = 0.0
        Nxyz = np.prod(self.domain.dims)
        while t < final_time:
            error_p = 1e5
            error_sw = 1e5
            loop_count = 0
            while error_p > self.eps_p or error_sw > self.eps_sw:
                loop_count += 1
                if loop_count > 10:
                    break
                pgrad = gradientTerm(self.pressure)
                sw_face = upwindMean(self.saturation, -pgrad)
                labdao = self.calculate_oil_mobility(sw_face)
                labdaw = self.calculate_water_mobility(sw_face)
                dlabdaodsw = self.calculate_oil_mobility_derivative(sw_face)
                dlabdawdsw = self.calculate_water_mobility_derivative(sw_face)
                labda = labdao + labdaw
                dlabdadsw = dlabdaodsw + dlabdawdsw
                Mdiffp1 = diffusionTerm(-labda)
                Mdiffp2 = diffusionTerm(-labdaw)
                Mconvsw1 = convectionUpwindTerm(-dlabdadsw * pgrad)
                Mconvsw2 = convectionUpwindTerm(-dlabdawdsw * pgrad)
                [Mtranssw2, RHStrans2] = transientTerm(self.initial_sw, dt, self.poros_field)
                RHS1 = divergenceTerm(-dlabdadsw * sw_face * pgrad)
                RHS2 = divergenceTerm(-dlabdawdsw * sw_face * pgrad)
                [Mbcp, RHSbcp] = boundaryConditionTerm(self.pressure_bc)
                [Mbcsw, RHSbcsw] = boundaryConditionTerm(self.saturation_bc)
                M = vstack([hstack([Mdiffp1 + Mbcp, Mconvsw1]), hstack([Mdiffp2, Mconvsw2 + Mtranssw2 + Mbcsw])])
                RHS = vstack([RHS1 + RHSbcp, RHS2 + RHStrans2 + RHSbcsw])
                x = spsolve(M, RHS)
                p_new = np.reshape(x[0:(Nxyz + 2)])
                sw_new = np.reshape(x[(Nxyz + 2):])
                error_p = max(abs((p_new - self.pressure.value[:]) / p_new))
                error_sw = max(abs(sw_new - self.saturation.value[:]))
                self.pressure.value[:] = p_new
                self.saturation.value[:] = sw_new
            if loop_count > 10:
                self.pressure.value[:] = self.initial_pressure.value[:]
                self.saturation.value[:] = self.initial_sw.value[:]
                dt = dt/5
                continue
            
            dsw = np.max(abs(sw_new[:]-self.initial_sw.value[:])/sw_new[:])
            t+=dt
            dt = np.min([dt*(self.dsw_allowed/dsw), 2*dt, final_time-t])
            self.initial_pressure.value[:] = self.pressure.value[:]
            self.initial_sw.value[:] = self.saturation.value[:]

if __name__ == "__main__":
    # define relative permeability parametrs
    krw0 = 1.0
    kro0 = 0.76
    nw = 2.4
    no = 2.0
    sor = 0.12
    swc = 0.09
    rel_perm = RelativePermeability(krw0=krw0, kro0=kro0, nw=nw, no=no, swc=swc, sor=sor)

    # define fluid properties
    mu_oil = 2e-3  # [Pa.s] oil viscosity
    mu_water = 1e-3  # [Pa.s] water viscosity
    phases = Fluids(mu_oil=mu_oil, mu_water=mu_water)

    # define reservoir properties
    k0 = 2e-12  # [m^2] average reservoir permeability
    phi0 = 0.2  # average porosity
    p0 = 100e5  # [bar] pressure
    sw0 = swc+0.1  # initial water saturation
    sw_in = 1
    field = Resevoir(rel_perm=rel_perm, fluids=phases, porosity=phi0, permeability=k0, 
                    sw_init=sw0, pressure_init=p0)

    # define operational conditions
    pin = 150e5  # [bar] injection pressure at the left boundary
    u_in = 1.0/(24*3600)  # [m/s] equal to 1 m/day
    p_back = p0  # [bar] production (back) pressure at the right boundary
    op_cond = OperationalConditions(injection_velocity=u_in, 
                                    injection_pressure=pin, production_pressure=p_back)

    # define the geometry
    Nx = 100  # number of cells in x direction
    W = 300  # [m] length of the domain in x direction
    m = createMesh1D(Nx, W)  # creates a 1D mesh

    # reservoir
    core = CoreModel1D(field=field, mesh=m, operational_conditions=op_cond)
    core.simulate(final_time=2*3600., dt=100.)
    # Lo = @(sw)(k/mu_oil*kro(sw))
    # dLwdsw = @(sw)(k/mu_water*dkrwdsw(sw))
    # dLodsw = @(sw)(k/mu_oil*dkrodsw(sw))
    # # Define the boundaries
    # BCp = createBC(m)  # Neumann BC for pressure
    # BCs = createBC(m)  # Neumann BC for saturation
    # # left boundary pressure gradient
    # BCp.left.a[:] = (krw(sw_in)*lw.xvalue(1, :)+kro(sw_in)*lo.xvalue(1, : )) BCp.left.b(: ) = 0 BCp.left.c(: ) = -u_in
    # # change the right boandary to constant pressure (Dirichlet)
    # # BCp.left.a[:]=0 BCp.left.b[:]=1 BCp.left.c[:]=pin
    # BCp.right.a[:] = 0 BCp.right.b(: ) = 1 BCp.right.c(: ) = p0
    # # change the left boundary to constant saturation (Dirichlet)
    # BCs.left.a[:] = 0 BCs.left.b(: ) = 1 BCs.left.c[:] = 1
    # # define the time step and solver properties
    # # dt = 1000 # [s] time step
    # dt = (W/Nx)/u_in/10  # [s]
    # t_end = 1000*dt  # [s] final time
    # eps_p = 1e-5  # pressure accuracy
    # eps_sw = 1e-5  # saturation accuracy
    # # define the variables
    # sw_old = createCellVariable(m, sw0, BCs)
    # p_old = createCellVariable(m, p0, BCp)
    # sw = sw_old
    # p = p_old
    # uw = -gradientTerm(p_old)  # an estimation of the water velocity
    # # start the main loop
    # # generate intial pressure profile (necessary to initialize the fully
    # # implicit solver)
    # dp_alwd = 100.0  # Pa
    # dsw_alwd = 0.05
    # t = 0
    # # fprintf(1, 'progress (##):  ')
    # while (t < t_end)
    # error_p = 1e5
    # error_sw = 1e5
    # # Implicit loop
    # loop_count = 0
    #     while ((error_p > eps_p) | | (error_sw > eps_sw))
    #     loop_count = loop_count+1
    #     if loop_count > 10
    #         break
    #         end
    #         # calculate parameters
    #         pgrad = gradientTerm(p)
    #         sw_face = upwindMean(sw, -pgrad)  # average value of water saturation
    #         labdao = lo.*funceval(kro, sw_face)
    #         labdaw = lw.*funceval(krw, sw_face)
    #         dlabdaodsw = lo.*funceval(dkrodsw, sw_face)
    #         dlabdawdsw = lw.*funceval(dkrwdsw, sw_face)
    #         labda = labdao+labdaw
    #         dlabdadsw = dlabdaodsw+dlabdawdsw
    #         # compute [Jacobian] matrices
    #         Mdiffp1 = diffusionTerm(-labda)
    #         Mdiffp2 = diffusionTerm(-labdaw)
    #         Mconvsw1 = convectionUpwindTerm(-dlabdadsw.*pgrad)
    #         Mconvsw2 = convectionUpwindTerm(-dlabdawdsw.*pgrad)
    #         [Mtranssw2, RHStrans2] = transientTerm(sw_old, dt, phi)
    #         # Compute RHS values
    #         RHS1 = divergenceTerm(-dlabdadsw.*sw_face.*pgrad)
    #         RHS2 = divergenceTerm(-dlabdawdsw.*sw_face.*pgrad)
    #         # include boundary conditions
    #         [Mbcp, RHSbcp] = boundaryCondition(BCp)
    #         [Mbcsw, RHSbcsw] = boundaryCondition(BCs)
    #         # Couple the equations BC goes into the block on the main diagonal
    #         M = [Mdiffp1+Mbcp Mconvsw1 Mdiffp2 Mconvsw2+Mtranssw2+Mbcsw]
    #         RHS = [RHS1+RHSbcp RHS2+RHStrans2+RHSbcsw]
    #         # solve the linear system of equations
    #         x = M\RHS
    #         # x = agmg(M, RHS, [], 1e-10, 500, [], [p.value[:] sw.value[:]])
    #         # separate the variables from the solution
    #         p_new = reshapeCell(m, full(x(1: (Nx+2))))
    #         sw_new = reshapeCell(m, full(x((Nx+2)+1: end)))
    #         # calculate error values
    #         error_p = max(abs((p_new(: )-p.value(: ))./p_new(: )))
    #         error_sw = max(abs(sw_new[:]-sw.value(: )))
    #         # assign new values of p and sw
    #         p.value = p_new
    #         sw.value = sw_new
    #     end
    #     if loop_count > 10
    #     p = p_old
    #     sw = sw_old
    #         dt = dt/5
    #         continue
    #     end
    #     dsw = max(abs(sw_new(: )-sw_old.value[:])./sw_new[:])
    #     t = t+dt
    #     fprintf(1, '\b\b\b\b\b\b\b\b\b\b\b\b\b\bProgress: #d ##', floor(t/t_end*100))
    #     dt = min([dt*(dsw_alwd/dsw), 2*dt, t_end-t])
    #     p_old = p
    #     sw_old = sw
    #     figure(1)visualizeCells(sw) drawnow
    # end
    # fprintf('\n')
        