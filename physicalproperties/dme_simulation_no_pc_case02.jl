using PyPlot, Polynomials, CoolProp, Roots, JLD
using JFVM
include("rel_perms_real.jl")

# some data
MW_water = 0.018 # [kg/mol]
MW_oil   = 0.142 # [kg/mol]
MW_DME   = 0.04607 # [kg/mol]

# Fig 4, SPE-177919, T = 70 degC, p = 2000 psi
x_DME_oil = [0.0, 0.04, 0.10, 0.39, 0.62]         # [mol frac]
w_DME_oil = MW_DME*x_DME_oil./(MW_DME*x_DME_oil+MW_oil*(1-x_DME_oil))
ρ_oil     = [0.748, 0.746, 0.744, 0.722, 0.697]*1000   # [g/cm^3]*1000 = [kg/m^3]

rho_fit = polyfit(w_DME_oil, ρ_oil, 2)

# plot(w_DME_oil, ρ_oil, "o", w_DME_oil, polyval(rho_fit, w_DME_oil))
# xlabel("DME mass fraction in oil [-]")
# ylabel("Oil density [g/cm^3]")

# experimental data from SPE-179771
# data for 1.4% brine and 71 degC and 2000 psia
x_DME_oil = [0.15, 0.30, 0.5, 0.6]   # [mol frac] [0.5, 13] is a madeup point
w_DME_oil = MW_DME*x_DME_oil./(MW_DME*x_DME_oil+MW_oil*(1-x_DME_oil)) # [mass frac]
K_value   = [17.0, 15.0, 13.0, 12.5]  # [-]
x_DME_water = x_DME_oil./K_value # [mole frac]
w_DME_water = MW_DME*x_DME_water./(MW_DME*x_DME_water+MW_water*(1-x_DME_oil))
K_value_mass = w_DME_oil./w_DME_water
K_fit = polyfit(w_DME_water, K_value_mass, 2)
w_DME_plot = collect(linspace(minimum(w_DME_water), maximum(w_DME_water), 50))
# plot(w_DME_water, K_value_mass, "o", w_DME_plot, polyval(K_fit,w_DME_plot))
# xlabel("DME mass fraction in water [-]")
# ylabel("K-value [-]")

# Fig 4, SPE-177919, T = 70 degC, p = 2000 psi
x_DME_oil = [0.0, 0.04, 0.10, 0.39, 0.62] # [mol frac]
w_DME_oil = MW_DME*x_DME_oil./(MW_DME*x_DME_oil+MW_oil*(1-x_DME_oil))
μ_oil = [1.04, 0.96, 0.84, 0.52, 0.34]/1000.0    # [cP]/1000 = [Pa.s]
mu_oil_fit = polyfit(w_DME_oil, μ_oil, 2)
w_DME_plot = collect(linspace(minimum(w_DME_oil), maximum(w_DME_oil), 50))
# plot(w_DME_oil, μ_oil, "o", w_DME_plot, polyval(mu_oil_fit, w_DME_plot))
# xlabel("DME mass fraction in oil [-]")
# ylabel("Oil viscosity [Pa.s]")

T0 = 70 + 273.15 # [K]
p0 = 2000/14.7*1e5   # [Pa]
mu_water = PropsSI("viscosity", "T", T0, "P", p0, "water") # [Pa.s]
mu_water_DME = [mu_water, 2*mu_water]
w_DME_water = [0.0, maximum(w_DME_oil)/maximum(K_value_mass)]
mu_water_fit = polyfit(w_DME_water, mu_water_DME, 1)
# plot(w_DME_water, mu_water_DME, "o", w_DME_water, mu_water_fit(w_DME_water))
# xlabel("DME mass fraction in water [-]")
# ylabel("Water-DME viscosity [Pa.s]")

# T = 70 + 273.15 # [K]
# p = 2000/14.7*1e5   # [Pa]
rho_water = PropsSI("D", "T", T0, "P", p0, "water")
# Corey rel-perm parameters
krw0 = 0.2
kro0 = 0.8
n_o  = 2.0
n_w  = 2.0
swc  = 0.08
sor  = 0.3
sorting_factor = 2.4
pce = 100 # [Pa]
pc0 = 1e5 # [Pa]
contact_angle = deg2rad(20) # [radian]

perm_val  = 0.01e-12 # [m^2] permeability
poros_val = 0.40     # [-] porosity

KRW  = sw -> krw.(sw, krw0, sor, swc, n_w)
dKRW = sw -> dkrwdsw.(sw, krw0, sor, swc, n_w)
KRO  = sw -> kro.(sw, kro0, sor, swc, n_o)
dKRO = sw -> dkrodsw.(sw, kro0, sor, swc, n_o)
# PC   = sw -> pc_imb(sw, pce, swc, sor, teta=contact_angle,
#                 labda=sorting_factor, b = 0.0, pc01=pc0, pc02=pc0)
# dPC  = sw -> dpc_imb(sw, pce, swc, sor, teta=contact_angle,
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
ρ_oil = w_DME_oil -> polyval(rho_fit, w_DME_oil)
dρ_oil = w_DME_oil -> polyval(polyder(rho_fit), w_DME_oil)
μ_water = w_DME_water -> polyval(mu_water_fit, w_DME_water)
dμ_water = w_DME_water -> polyval(polyder(mu_water_fit), w_DME_water)
μ_oil = w_DME_oil -> polyval(mu_oil_fit, w_DME_oil)
dμ_oil = w_DME_oil -> polyval(polyder(mu_oil_fit), w_DME_oil)
Kval = w_DME_water -> polyval(K_fit, w_DME_water)
dKval = w_DME_water -> polyval(polyder(K_fit), w_DME_water)
w_oil = w_DME_water -> polyval(Poly([0.0,1.0])*K_fit, w_DME_water)
dw_oil = w_DME_water -> polyval(polyder(Poly([0.0,1.0])*K_fit), w_DME_water)

sw_plot = collect(linspace(0,1,1000))
# plot(sw_plot, KRW.(sw_plot), sw_plot, KRO.(sw_plot))
# axis([0,1,0,1])

# plot(sw_plot, PC.(sw_plot))
# plot([swc, swc], [-pc0, pc0], [1-sor, 1-sor], [-pc0, pc0],
#     [0,1], [0,0], "--")
# xlabel("Water saturation (Sw) [-]")
# ylabel("Capillary pressure (Pc) [Pa]")
sw_imb_end = fzero(PC, [swc, 1-sor])
# plot(sw_imb_end, 0, "o")
# axis([0, 1, -pc0, pc0])
# yscale("symlog")

Lx   = 0.3 # [m] Length of the core
Nx  = 50  # number of grids in the x direction
m   = createMesh1D(Nx, Lx)

u_inj = 1.0/(3600*24) # 1 m/day to m/s injection velocity
c_inj1 = 0.0 # DME mass fraction in the injected water before water flood
c_inj2 = 0.2 # DME mass fraction in the injected water after water flood

BCp = createBC(m) # pressure boundary condition
BCs = createBC(m) # saturation boundary condition
BCc = createBC(m) # concentration (DME in water) boundary condition

BCp.right.a[:] = 0.0
BCp.right.b[:] = 1.0
BCp.right.c[:] = p0

BCp.left.a[:]  = perm_val/μ_water(c_inj1)
BCp.left.b[:]  = 0.0
BCp.left.c[:]  = -u_inj

BCs.left.a[:]  = 0.0
BCs.left.b[:]  = 1.0
BCs.left.c[:]  = 1.0

# BCs.right.a[:]  = 0.0
# BCs.right.b[:]  = 1.0
# BCs.right.c[:]  = sw_imb_end

BCc.left.a[:]  = 0.0
BCc.left.b[:]  = 1.0
BCc.left.c[:]  = c_inj1

# discretize
M_bc_p, RHS_bc_p = boundaryConditionTerm(BCp)
M_bc_s, RHS_bc_s = boundaryConditionTerm(BCs)
M_bc_c, RHS_bc_c = boundaryConditionTerm(BCc)

c0  = 0.0 # [mass frac] initial DME concentration in the water phase
sw0 = swc+0.05 # [vol frac] initial water saturation

p_init  = createCellVariable(m, p0, BCp)
c_init  = createCellVariable(m, c0, BCc)
sw_init = createCellVariable(m, sw0, BCs)

# new values of each variable
p_val  = createCellVariable(m, p0, BCp)
c_val  = createCellVariable(m, c0, BCc)
sw_val = createCellVariable(m, sw0, BCs)

Δsw_init = sw_init - sw_val # we solve for this variable
Δc_init  = c_init - c_val   # we solve for this variable

# Other cell variables
k = createCellVariable(m, perm_val)
ϕ = createCellVariable(m, poros_val)

n_pv1     = 0.1 # number of injected pore volumes
t_final1  = n_pv1*Lx/(u_inj/poros_val) # [s] final time for water flood
dt01      = t_final1/n_pv1/Nx # [s] time step
dt        = dt01
n_pv2     = 0.1 # with DME
t_final2  = t_final1+n_pv1*Lx/(u_inj/poros_val) # [s] final time for DME flood
dt02      = t_final2/n_pv2/Nx/10 # [s] time step
dt0 = dt01

# outside the loop: initialization
uw       = gradientTerm(p_val) # only for initialization of the water velocity vector
k_face   = harmonicMean(k)     # permeability on the cell faces

# this whole thing goes inside two loops (Newton and time)
tol_s = 1e-7
tol_c = 1e-7
max_change_c = 0.1 # 10 % relative change
max_change_s = 0.1 # 10 % relative change
max_int_loop = 4
t = dt
oil_init=domainInt(1-sw_init) # initial oil volume in the core
rec_fact=zeros(1)
dp_hist = zeros(1)
t_s=zeros(1)

while t<t_final2
  if t>t_final1
    dt0 = dt02
    BCc.left.a[:]  = 0.0
    BCc.left.b[:]  = 1.0
    BCc.left.c[:]  = c_inj2
    BCp.left.a[:]  = perm_val/μ_water(c_inj2)
    BCp.left.b[:]  = 0.0
    BCp.left.c[:]  = -u_inj
    c_init  = createCellVariable(m, internalCells(c_val), BCc)
    p_val   = createCellVariable(m, internalCells(p_val), BCp)
    M_bc_p, RHS_bc_p = boundaryConditionTerm(BCp)
    M_bc_c, RHS_bc_c = boundaryConditionTerm(BCc)
  end
  error_s = 2*tol_s
  error_c = 2*tol_c
  loop_countor = 0
  while error_s>tol_s || error_c>tol_c
    loop_countor += 1
    if loop_countor > max_int_loop
      sw_val = copyCell(sw_init)
      c_val  = copyCell(c_init)
      p_val  = copyCell(p_init)
      dt = dt/3.0
      break
    end
    ∇p0      = gradientTerm(p_val)
    ∇s0      = gradientTerm(sw_val)

    c_face   = upwindMean(c_val, uw)
    c_oil         = cellEval(w_oil, c_val) # [mass frac] DME concentration in oil
    c_oil_face    = faceEval(w_oil, c_face)
    sw_face       = upwindMean(sw_val, uw)

    mu_water_face = faceEval(μ_water, c_face)
    dμ_water_face = faceEval(dμ_water, c_face)
    mu_oil_face   = faceEval(μ_oil, c_oil_face)
    dμ_oil_face   = faceEval(dw_oil, c_face).*faceEval(dμ_oil,c_oil_face)
    rho_oil_cell  = cellEval(ρ_oil, c_oil)
    rho_oil_face  = faceEval(ρ_oil, c_oil_face)
    dρ_oil_cell   = cellEval(dw_oil, c_val).*cellEval(dρ_oil, c_oil)
    dρ_oil_face   = faceEval(dw_oil, c_face).*faceEval(dρ_oil, c_oil_face)
    dρc_oil_cell = dρ_oil_cell.*c_oil+cellEval(dw_oil, c_val).*rho_oil_cell
    dρc_oil_face = dρ_oil_face.*c_oil_face+faceEval(dw_oil, c_face).*rho_oil_face
    dpc_face      = faceEval(dPC, sw_face)
    d2pc_face      = faceEval(d2PC, sw_face)

    dρ_dμ_oil = (dρ_oil_face.*mu_oil_face-dμ_oil_face.*rho_oil_face)./(mu_oil_face.*mu_oil_face)
    dρc_dμ_oil = (dρc_oil_face.*mu_oil_face-dμ_oil_face.*rho_oil_face.*c_oil_face)./(mu_oil_face.*mu_oil_face)

    krw_face      = faceEval(KRW,sw_face)
    kro_face      = faceEval(KRO,sw_face)

    λ_w_face      = k_face.*krw_face./mu_water_face
    λ_o_face      = k_face.*kro_face./mu_oil_face
    uw            = -λ_w_face.*∇p0

    dλ_w_face     = k_face.*faceEval(dKRW,sw_face)./mu_water_face
    dλ_o_face     = k_face.*faceEval(dKRO,sw_face)./mu_oil_face

    # λ_w_face     = harmonicMean(λ_w) # Wrong!
    # λ_o_face = harmonicMean(λ_o)     # Wrong!

    # Δsw_init = sw_init - sw_val # we solve for this variable
    # Δc_init  = c_init - c_val   # we solve for this variable

    # water mass balance (wmb)
    M_t_s_wmb, RHS_t_s_wmb = transientTerm(sw_init, dt, rho_water*ϕ)
    M_d_p_wmb = diffusionTerm(-rho_water*λ_w_face)
    M_a_s_wmb = convectionUpwindTerm(-rho_water*dλ_w_face.*∇p0)
    M_a_c_wmb = convectionUpwindTerm(rho_water*dμ_water_face./mu_water_face.*λ_w_face.*∇p0)

    # oil mass balance (omb)
    M_t_s_omb, RHS_t_s_omb = transientTerm(sw_init, dt, -rho_oil_cell.*ϕ)
    M_t_c_omb, RHS_t_c_omb = transientTerm(c_init, dt, (1.0-sw_val).*dρ_oil_cell.*ϕ)
    M_d_p_omb = diffusionTerm(-rho_oil_face.*λ_o_face)
    # M_a_c_omb = convectionUpwindTerm(-dρ_dμ_oil.*k_face.*kro_face.*(∇p0+dpc_face.*∇s0))
    M_a_c_omb = convectionUpwindTerm(-dρ_dμ_oil.*k_face.*kro_face.*∇p0)
    M_a_s_omb = convectionUpwindTerm(-rho_oil_face.*(dλ_o_face.*∇p0))
                                    # +(dλ_o_face.*dpc_face+λ_o_face.*d2pc_face).*∇s0))
    M_d_s_omb = 0.0 #diffusionTerm(-rho_oil_face.*λ_o_face.*dpc_face)

    # DME mass balance (dme)
    M_t_s_dme, RHS_t_s_dme = transientTerm(sw_init, dt, (rho_water*c_val-rho_oil_cell.*c_oil).*ϕ)
    M_t_c_dme, RHS_t_c_dme = transientTerm(c_init, dt, (rho_water*sw_val+
                                           (1.0-sw_val).*dρc_oil_cell).*ϕ)
    M_d_p_dme  = diffusionTerm(-rho_water*c_face.*λ_w_face-rho_oil_face.*c_oil_face.*λ_o_face)
    M_a_s_dme  = convectionUpwindTerm(-(rho_water*c_face.*dλ_w_face+
                                      rho_oil_face.*c_oil_face.*dλ_o_face).*∇p0)
                                    #   -rho_oil_face.*c_oil_face.*
                                    #   (dλ_o_face.*dpc_face+λ_o_face.*d2pc_face).*∇s0)
    M_a_c_dme  = convectionUpwindTerm(-(rho_water*(mu_water_face-dμ_water_face.*c_face)
                                      ./mu_water_face.*λ_w_face+
                                        dρc_dμ_oil.*k_face.*kro_face).*∇p0)
                                    #   -dρc_dμ_oil.*k_face.*kro_face.*dpc_face.*∇s0)
    M_d_s_dme  = 0.0 # diffusionTerm(-c_oil_face.*rho_oil_face.*λ_o_face.*dpc_face)

    # create the PDE system M [p;s;c]=RHS
    # x_val = [p_val.value[:]; sw_val.value[:]; c_val.value[:]]
    M = [M_bc_p+M_d_p_wmb   M_t_s_wmb+M_a_s_wmb    M_a_c_wmb;
         M_d_p_omb   M_bc_s+M_t_s_omb+M_a_s_omb+M_d_s_omb    M_t_c_omb+M_a_c_omb;
         M_d_p_dme    M_t_s_dme+M_a_s_dme+M_d_s_dme    M_bc_c+M_t_c_dme+M_a_c_dme;]

    RHS = [RHS_bc_p+RHS_t_s_wmb+M_a_s_wmb*sw_val.value[:]+M_a_c_wmb*c_val.value[:];
           RHS_bc_s+RHS_t_s_omb+RHS_t_c_omb+(M_a_s_omb+M_d_s_omb)*sw_val.value[:]+M_a_c_omb*c_val.value[:];
           RHS_bc_c+RHS_t_s_dme+RHS_t_c_dme+(M_a_s_dme+M_d_s_dme)*sw_val.value[:]+M_a_c_dme*c_val.value[:];]

    # x_sol = solveLinearPDE(m, M, RHS)
    x_sol = M\RHS

    p_new = x_sol[1:Nx+2]
    s_new = x_sol[Nx+3:2*Nx+4]
    c_new = x_sol[2*Nx+5:end]
    error_s = sumabs(s_new[2:end-1]-sw_val.value[2:end-1])
    error_c = sumabs(c_new[2:end-1]-c_val.value[2:end-1])
    # println(error_s)
    # println(error_c)
    p_val.value[:] = p_new[:]
    sw_val.value[:] = s_new[:]
    c_val.value[:] = c_new[:]
  end
  if loop_countor<max_int_loop
    p_init = copyCell(p_val)
    sw_init = copyCell(sw_val)
    c_init = copyCell(c_val)
    t +=dt
    dt = dt0
    println("progress is $(t/t_final2*100) [%]")
    rec_fact=push!(rec_fact, (oil_init-domainInt(1-sw_val))/oil_init)
    t_s=push!(t_s, t)
    dp_hist=push!(dp_hist, 0.5*sum(p_val.value[1:2])-p0)
  end
end
figure(1)
visualizeCells(sw_init)
visualizeCells(c_init)
savefig("profile.png")
figure(2)
plot(t_s, rec_fact)
savefig("recovery.png")
figure(3)
plot(t_s, dp_hist)
savefig("dp.png")

# @save "filename.jld"
