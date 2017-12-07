using PyPlot, Polynomials, CoolProp, Roots
using JFVM
import GR
include("../functions/rel_perms_real.jl")

T0 = 70 + 273.15 # [K]
p0 = 2000/14.7*1e5   # [Pa]
mu_water = PropsSI("viscosity", "T", T0, "P", p0, "water") # [Pa.s]
rho_water = PropsSI("D", "T", T0, "P", p0, "water")
mu_oil    = 5e-3 # [Pa.s] oil viscosity
rho_oil   = 800  # [kg/m3] oil density
# Corey rel-perm parameters
krw0 = 0.2
kro0 = 0.8
n_o  = 2.0
n_w  = 2.0
swc  = 0.08
sor  = 0.3
sorting_factor = 2.4
pce = 1e5 # [Pa]
pc0 = 50*pce # [Pa]
contact_angle = deg2rad(20) # [radian]

perm_val  = 0.01e-12 # [m^2] permeability
poros_val = 0.40     # [-] porosity

KRW  = sw -> krw.(sw, krw0, sor, swc, n_w)
dKRW = sw -> dkrwdsw.(sw, krw0, sor, swc, n_w)
KRO  = sw -> kro.(sw, kro0, sor, swc, n_o)
dKRO = sw -> dkrodsw.(sw, kro0, sor, swc, n_o)
# PC   = sw -> pc_imb.(sw, pce, swc, sor, teta=contact_angle,
#                 labda=sorting_factor, b = 0.0, pc01=pc0, pc02=pc0)
# dPC  = sw -> dpc_imb.(sw, pce, swc, sor, teta=contact_angle,
#                 labda1=sorting_factor, labda2=sorting_factor,
#                 b = 1.0, pc01=pc0, pc02=pc0)
# d2PC = sw -> d2pc_imb.(sw, pce, swc, sor, teta=contact_angle,
#                 labda1=sorting_factor, labda2=sorting_factor,
#                 b = 1.0, pc01=pc0, pc02=pc0)
# PC2  = sw -> pc_imb2(sw, pce, swc, sor, teta=contact_angle,
#                 labda=sorting_factor, b = 0.0, pc_star1=pc0,
#                 pc_star2=pc0)
# dPC2 = sw -> dpc_imb2(sw, pce, swc, sor, teta=contact_angle,
#                 labda1=sorting_factor, labda2=sorting_factor,
#                 b = 1.0, pc_star1=pc0, pc_star2=pc0)
PC  = sw -> pc_imb3.(sw, pce, swc, sor, teta=contact_angle,
                labda=sorting_factor, b = 0.0, pc_star1=pc0,
                pc_star2=pc0)
dPC = sw -> dpc_imb3.(sw, pce, swc, sor, teta=contact_angle,
                labda1=sorting_factor, labda2=sorting_factor,
                b = 1.0, pc_star1=pc0, pc_star2=pc0)
d2PC = sw -> d2pc_imb3.(sw, pce, swc, sor, teta=contact_angle,
                labda1=sorting_factor, labda2=sorting_factor,
                b = 1.0, pc_star1=pc0, pc_star2=pc0)
# PCdrain2 = sw -> pc_drain2(sw, pce, swc, labda=sorting_factor, pc_star=pc0)
# PCdrain3 = sw -> pc_drain3(sw, pce, swc, labda=sorting_factor, pc_star=pc0)
sw_imb_end = fzero(PC, [swc, 1-sor])

Lx   = 0.3 # [m]
Nx  = 500  # number of grids in the x direction
m   = createMesh1D(Nx, Lx)

u_inj = 1.0/(3600*24) # 1 m/day to m/s injection velocity

BCp = createBC(m) # pressure boundary condition
BCs = createBC(m) # saturation boundary condition

BCp.right.a[:] = 0.0
BCp.right.b[:] = 1.0
BCp.right.c[:] = p0

BCp.left.a[:]  = perm_val/mu_water
BCp.left.b[:]  = 0.0
BCp.left.c[:]  = -u_inj

BCs.left.a[:]  = 0.0
BCs.left.b[:]  = 1.0
BCs.left.c[:]  = 1.0

BCs.right.a[:]  = 0.0
BCs.right.b[:]  = 1.0
BCs.right.c[:]  = sw_imb_end
boundary_switch = true

# discretize
M_bc_p, RHS_bc_p = boundaryConditionTerm(BCp)
M_bc_s, RHS_bc_s = boundaryConditionTerm(BCs)

sw0 = swc+0.15 # [vol frac] initial water saturation

p_init  = createCellVariable(m, p0, BCp)
sw_init = createCellVariable(m, sw0, BCs)

# new values of each variable
p_val  = createCellVariable(m, p0, BCp)
sw_val = createCellVariable(m, sw0, BCs)

Δsw_init = sw_init - sw_val # we solve for this variable

# Other cell variables
k = createCellVariable(m, perm_val)
ϕ = createCellVariable(m, poros_val)

n_pv    = 0.3 # number of injected pore volumes
t_final = n_pv*Lx/(u_inj/poros_val) # [s] final time
dt0      = t_final/n_pv/Nx/100 # [s] time step
dt = dt0
# outside the loop: initialization
ut       = -gradientTerm(p_val) # only for initialization of the water velocity vector
k_face   = harmonicMean(k)     # permeability on the cell faces

# this whole thing goes inside two loops (Newton and time)
tol_s = 1e-7
max_change_s = 0.1 # 10 % relative change
max_int_loop = 10
t = dt
while t<t_final
    error_s = 2*tol_s
    loop_countor = 0
    while error_s>tol_s
        loop_countor += 1
        if loop_countor > max_int_loop
          sw_val = copyCell(sw_init)
          p_val  = copyCell(p_init)
          dt = dt/10.0
          break
        end
        ∇p0      = gradientTerm(p_val)
        ∇s0      = gradientTerm(sw_val)

        sw_face       = upwindMean(sw_val, ut)

        dpc_face      = faceEval(dPC, sw_face)
        d2pc_face      = faceEval(d2PC, sw_face)

        krw_face      = faceEval(KRW,sw_face)
        kro_face      = faceEval(KRO,sw_face)

        λ_w_face      = k_face.*krw_face/mu_water
        λ_o_face      = k_face.*kro_face/mu_oil
        uw            = -λ_w_face.*∇p0 #-λ_w_face.*(∇p0+dpc_face.*∇s0)
        uo            = -λ_w_face.*(∇p0+dpc_face.*∇s0)
        ut = uw+uo

        dλ_w_face     = k_face.*faceEval(dKRW,sw_face)/mu_water
        dλ_o_face     = k_face.*faceEval(dKRO,sw_face)/mu_oil

        # λ_w_face     = harmonicMean(λ_w) # Wrong!
        # λ_o_face = harmonicMean(λ_o)     # Wrong!

        # Δsw_init = sw_init - sw_val # we solve for this variable
        # Δc_init  = c_init - c_val   # we solve for this variable

        # water mass balance (wmb)
        M_t_s_wmb, RHS_t_s_wmb = transientTerm(sw_init, dt, rho_water*ϕ)
        M_d_p_wmb = diffusionTerm(-rho_water*λ_w_face)
        M_a_s_wmb = convectionTerm(-rho_water*dλ_w_face.*∇p0)

        # oil mass balance (omb)
        M_t_s_omb, RHS_t_s_omb = transientTerm(sw_init, dt, -rho_oil*ϕ)
        M_d_p_omb = diffusionTerm(-rho_oil*λ_o_face)
        M_a_s_omb = convectionTerm(-rho_oil*dλ_o_face.*(∇p0+dpc_face.*∇s0))
        RHS_d_s_omb = divergenceTerm(rho_oil*λ_o_face.*dpc_face.*∇s0)
        # M_a_s_omb = convectionUpwindTerm(-rho_oil*((dλ_o_face.*∇p0)
        #                                 +(dλ_o_face.*dpc_face+λ_o_face.*d2pc_face).*∇s0), uw)
        M_d_s_omb = 0.0 #diffusionTerm(-rho_oil*λ_o_face.*dpc_face)
        
        # create the PDE system M [p;s;c]=RHS
        # x_val = [p_val.value[:]; sw_val.value[:]; c_val.value[:]]
        M = [M_bc_p+M_d_p_wmb   M_t_s_wmb+M_a_s_wmb;
             M_d_p_omb   M_bc_s+M_t_s_omb+M_a_s_omb+M_d_s_omb;]

        RHS = [RHS_bc_p+RHS_t_s_wmb+M_a_s_wmb*sw_val.value[:];
               RHS_bc_s+RHS_t_s_omb+RHS_d_s_omb+(M_a_s_omb)*sw_val.value[:];]
        # x_sol = solveLinearPDE(m, M, RHS)
        x_sol = M\RHS

        p_new = x_sol[1:Nx+2]
        s_new = x_sol[Nx+3:2*Nx+4]
        error_s = sum(abs, s_new[2:end-1]-sw_val.value[2:end-1])
        # println(error_s)
        p_val.value[:] = p_new[:]
        dsw_apple = 0.2
        eps_apple = sqrt(eps()) # 1e-5
        for i in eachindex(s_new)
            if s_new[i]>=(1-sor)
                if sw_val.value[i]<(1-sor-eps_apple)
                    sw_val.value[i] = 1-sor-eps_apple
                else
                    sw_val.value[i] = 1-sor
                end
            elseif s_new[i]<=swc
                if sw_val.value[i]> swc+eps_apple
                    sw_val.value[i] = swc+eps_apple
                else
                    sw_val.value[i] = swc
                end
            elseif abs(s_new[i]-sw_val.value[i])>dsw_apple
                sw_val.value[i] += dsw_apple*sign(s_new[i]-sw_val.value[i])
            else
                sw_val.value[i] = s_new[i]
            end
        end
        sw_val = createCellVariable(m, internalCells(sw_val), BCs)
        # sw_val.value[:] = s_new[:]
    end
    if loop_countor<max_int_loop
      p_init = copyCell(p_val)
      sw_init = copyCell(sw_val)
      t +=dt
      dt = dt0
      GR.plot(sw_init.value[2:end-1])
      println("progress: $(t/t_final*100) %")
      if 0.5*sum(sw_val.value[end-1:end])>=sw_imb_end && boundary_switch
          boundary_switch = false
          BCs.right.a[:]  = 0.0
          BCs.right.b[:]  = 1.0
          BCs.right.c[:]  = sw_imb_end
          M_bc_s, RHS_bc_s = boundaryConditionTerm(BCs)
      end
    end
end
visualizeCells(sw_init)
# savefig("profile.png")
