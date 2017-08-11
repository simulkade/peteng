"""
Buckley Leverett equation
dependent variables: water saturation
Prepared for educational purposes by Ali A Eftekhari
"""
function forced_imb_impes(mu_water, mu_oil, ut, phi0,
  perm_val, swc, sor, kro0, no, krw0, nw, sw0, sw_inj, L, pv_inj; Nx = 100)
  m = createMesh1D(Nx, L)
  # m = createMeshCylindrical2D(Nx, Ny, W, H) # creates a 2D mesh
  ## define the physical parametrs

  p0 = 1e5 # [bar] pressure
  pin = 150e5 # [bar] injection pressure at the left boundary
  u_in= ut # [m/s] equal to 1 m/day
  sw_in = sw_inj

  eps1=1e-6

  k=createCellVariable(m, perm_val)
  phi=createCellVariable(m, phi0)

  lw = geometricMean(k)/mu_water
  lo = geometricMean(k)/mu_oil

  ## Define the boundaries: all fixed Sw=1, fixed pressure everywhere(?)
  BCp = createBC(m) # Neumann BC for pressure
  BCs = createBC(m) # Neumann BC for saturation
  # left boundary pressure gradient
  # change the right boandary to constant pressure (Dirichlet)
  BCp.left.a[:]=perm_val/mu_water
  BCp.left.b[:]=0
  BCp.left.c[:]=-u_in
  BCp.right.a[:]=0.0
  BCp.right.b[:]=1.0
  BCp.right.c[:]=p0

  BCs.left.a[:]=0.0
  BCs.left.b[:]=1.0
  BCs.left.c[:]=sw_in

  t_end = pv_inj*phi0*L/u_in # [s] final time
  dt=L/Nx/u_in
  eps_p = 1e-7 # pressure accuracy
  eps_sw = 1e-7 # saturation accuracy
  ## define the variables
  sw_old = createCellVariable(m, sw0, BCs)
  p_old = createCellVariable(m, p0, BCp)
  sw = copyCell(sw_old)
  oil_init=domainInt(1-sw_old)
  p = copyCell(p_old)
  uw = -gradientTerm(p_old) # an estimation of the water velocity
  ## start the main loop
  # generate intial pressure profile (necessary to initialize the fully
  # implicit solver)
  rec_fact=zeros(1)
  t_day=zeros(1)
  dp_hist = zeros(1)
  t = 0.0
  dt0=dt
  dsw_alwd= 0.05
  dp_alwd= 100.0 # Pa
  labdaw=createFaceVariable(m,1.0) # initialize
  # prog_1=ProgressThresh(t_end)
  while (t<t_end)
    print("progress is $((t/t_end*100)) %", "\u1b[1G")
    # for i=1:5
      error_p = 1e5
      error_sw = 1e5
      # Implicit loop
  #     while ((error_p>eps_p) || (error_sw>eps_sw))
      while(true)
          # calculate parameters
          pgrad = gradientTerm(p)
          sw_face = upwindMean(sw, -labdaw.*pgrad) # average value of water saturation
          sw_ave=arithmeticMean(sw)
          # solve for pressure at known Sw
          labdao = lo.*faceEval(s -> kro.(s, kro0, sor, swc, no), sw_face)
          labdaw = lw.*faceEval(s -> krw.(s, krw0, sor, swc, nw), sw_face)
          labda = labdao+labdaw
          # compute [Jacobian] matrices
          Mdiffp1 = diffusionTerm(-labda)
          Mbcp, RHSbcp = boundaryConditionTerm(BCp)
          RHS1 = RHSbcp # with capillary
          p_new=solvePDE(m, Mdiffp1+Mbcp, RHS1)

          # solve for Sw
          pgrad = gradientTerm(p_new)
          uw=-labdaw.*pgrad
          Mbcsw, RHSbcsw = boundaryConditionTerm(BCs)
          RHS_sw=-divergenceTerm(uw)
          sw_new=solveExplicitPDE(sw_old, dt, RHS_sw, BCs, phi)

          error_p = maximum(abs((internalCells(p_new)-internalCells(p))./internalCells(p_new)))
          error_sw = maximum(abs(internalCells(sw_new)-internalCells(sw)))
          dt_new=dt*min(dp_alwd/error_p, dsw_alwd/error_sw)

          # assign new values of p and sw
          if error_sw>dsw_alwd
              dt=dt*(dsw_alwd/error_sw)
              # print("new time step = $dt \n")
          else
              t=t+dt
              p = copyCell(p_new)
              sw = copyCell(sw_new)
              p_old = copyCell(p)
              sw_old = copyCell(sw)
              dt=min(dt*(dsw_alwd/error_sw), 15*dt, t_end-t-eps())
              break
          end
      end

      push!(rec_fact, (oil_init-domainInt(1-sw))/oil_init)
      push!(t_day, t)
      push!(dp_hist, 0.5*sum(p.value[1:2])-p0)
      # print(t)
      # GR.plot(sw.value[2:end-1])
      # GR.imshow(sw.value[2:end-1,2:end-1])
      # visualizeCells(1-sw)
      # GR.plot(t_day/3600/24, rec_fact)
      # xlabel('time [day]')
      # ylabel('recovery factor')
      # title([num2str(t/3600/24) ' day']) drawnow
  end
  dp_hist[1] = dp_hist[2]
  t_day, rec_fact, dp_hist, sw
end # imb_impes
