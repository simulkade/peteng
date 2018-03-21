using JFVM, PyPlot, JPhreeqc

# Physical properties
n_comp = 6 # number of components
Di = [1.0e-19, 2.0e-9, 3.0e-9] # m^2/s
u_inj = 1.0/(24*3600) # m/s
c_init = [0.0, 0.5, 0.1] # mol/m^3
phi_val = 0.4 # [-] porosity
k_val = 0.001e-12 # m^2 permeability
c_bound = [1.0, 0.0, 2.0] # mol/m^3 left boundary

# Define the domain
Lx = 1.0 # m
Nx = 100 # [-]
nxyz = Nx
m = createMesh1D(Nx, Lx)
v_cell = cellVolume(m)
sw = createCellVariable(m, 1.0)
phi = createCellVariable(m, phi_val)

# m = createMesh2D(Nx, 10, Lx, Lx/10)

nthreads=1 # number of threads
id=RM_Create(nxyz, nthreads)
setDefaultPhreeqcProperties(id) # set properties
setDefaultPhreeqcUnits(id)
status = RM_SetUnitsSurface(id, 0)       # 0, mol/L cell 1, mol/L water 2 mol/L rock
status = RM_SetUnitsSolution(id, 2)
status = RM_UseSolutionDensityVolume(id, 2)
status = RM_SetSelectedOutputOn(id, 1)       # enable selected output
grid2chem = collect(Int, 0:nxyz-1)
status = RM_CreateMapping(id, grid2chem)
status = RM_SetPrintChemistryOn(id, 0, 1, 0) # workers, initial_phreeqc, utility

RM_LoadDatabase(id, "PHREEQC.DAT")
vol=internalCells(v_cell)[:]*1000 # cel volume in liter
RM_SetRepresentativeVolume(id, vol)
por=internalCells(phi)[:]
RM_SetPorosity(id, por)
# Set initial saturation
sat = internalCells(sw)[:]
status = RM_SetSaturation(id, sat)
RM_SetDensity(id, ones(nxyz))
RM_RunFile(id, 1,1,1, "phreeqc_input.pqi")

RM_FindComponents(id)

ic1 = zeros(Int, nxyz, 7)
ic2 = zeros(Int, nxyz, 7)
f11 = zeros(Float64, nxyz, 7)
ic1[:,1] = 1:nxyz       # Solution 1
ic1[:,2] = -1      # Equilibrium phases none
ic1[:,3] = -1       # Exchange none
ic1[:,4] = -1      # Surface 1
ic1[:,5] = -1      # Gas phase none
ic1[:,6] = -1      # Solid solutions none
ic1[:,7] = -1      # Kinetics none
ic2[:,1] = -1      # Solution none
ic2[:,2] = -1      # Equilibrium phases none
ic2[:,3] = -1      # Exchange none
ic2[:,4] = -1      # Surface none
ic2[:,5] = -1      # Gas phase none
ic2[:,6] = -1      # Solid solutions none
ic2[:,7] = -1      # Kinetics none
f11[:,1]  = 1.0      # Mixing fraction ic1 Solution
f11[:,2] = 1.0      # Mixing fraction ic1 Equilibrium phases
f11[:,3] = 1.0      # Mixing fraction ic1 Exchange 1
f11[:,4] = 1.0      # Mixing fraction ic1 Surface
f11[:,5] = 1.0      # Mixing fraction ic1 Gas phase
f11[:,6] = 1.0      # Mixing fraction ic1 Solid solutions
f11[:,7] = 1.0      # Mixing fraction ic1 Kinetics
status = RM_InitialPhreeqc2Module(id, ic1[:], ic2[:], f11[:])

time1 = 0.0 # changed time to time1 to avoid conflicts with julia time function
time_step = 0.0
status = RM_SetTime(id, time1)
status = RM_SetTimeStep(id, time_step)
# get elements
# status = RM_SetSpeciesSaveOn(id, 0)
status = RM_RunCells(id)
element_list= getComponentList(id) # include charge
n_element= RM_GetComponentCount(id) # include charge
c = zeros(Float64, n_element*nxyz)
status = RM_GetConcentrations(id, c) # get the aqueous concentrations of elements
c=reshape(c, nxyz, n_element)

c_left=zeros(Float64, n_element)
n_bound=1
status= RM_InitialPhreeqc2Concentrations(id, c_left, n_bound, zeros(Int, n_bound),
  zeros(Int, n_bound), ones(Float64, n_bound)) # solution 0

# boundary conditions
bc = Array{BoundaryCondition}(n_comp)
for i in 1:n_comp
  bc[i] = createBC(m)
  bc[i].left.a[:] = 0.0
  bc[i].left.b[:] = 1.0
  bc[i].left.c[:] = c_left[i]
end

# Initial conditions
c_old = Array{CellValue}(n_comp)
c_new = Array{CellValue}(n_comp)
for i in 1:n_comp
  c_old[i] = createCellVariable(m, c[:,i][:], bc[i])
end

# velocity and diffusion coefficients
# D = Array{FaceValue}(n_comp)
u = createFaceVariable(m, u_inj)
D = createFaceVariable(m, Di[1])

# solver options
t_final = 0.5*Lx/(u_inj/phi_val)
dt = Lx/(u_inj/phi_val)/Nx

# Discretization
# M_diff = Array{Any}(3)
M_conv = convectionTerm(u)
M_diff = diffusionTerm(D)
M_bc = Array{Any}(n_comp)
RHS_bc = Array{Any}(n_comp)
for i in 1:n_comp
  M_bc[i], RHS_bc[i] = boundaryConditionTerm(bc[i])
end

M_t = Array{Any}(n_comp)
RHS_t = Array{Any}(n_comp)
c_reac= zeros(Float64, nxyz, n_element)

for t in dt:dt:t_final
  for i in 1:n_comp
    M_t[i], RHS_t[i] = transientTerm(c_old[i], dt, phi_val)
    c_new[i] = solveLinearPDE(m, M_t[i]+M_bc[i]+M_conv-phi_val*M_diff, RHS_bc[i]+RHS_t[i])
    c_reac[:,i]=internalCells(c_new[i])[:]
  end

  status= RM_SetConcentrations(id, c_reac[:])
  sat = internalCells(sw)[:]
  status = RM_SetSaturation(id, sat)
  status = RM_SetTime(id, t)
  status = RM_SetTimeStep(id, dt)
  status = RM_RunCells(id)
  c = zeros(Float64, n_element*nxyz)
  status = RM_GetConcentrations(id, c)

  c=reshape(c, nxyz, n_element)
  # update the concentrations
  # conc_init=Array{Any}(n_element+1)
  for i in 1:n_element
    c_old[i]=createCellVariable(m, c[:,i][:], bc[i])
  end
  println("progress: $(t/t_final*100) %")
end
figure(4)
visualizeCells(c_old[4])
visualizeCells(c_old[5])
visualizeCells(c_old[6])
savefig("profiles1.png")
