"""
Polymer flooding model
injecting polymer/water mixture into a reservoir
The adsorption of polymer and change of the injected water viscosity due to the 
change of concentration is included in the model.
Note: currently, the adsorption term is not included; it will be added later
"""
function polymer_implicit_upwind(core_props, fluids_pw, fluids_fw, rel_perms_fw,
  rel_perms_pw, core_flood; pv_inj = 0.3)
  c_insitu = 1.0
  c_inj = 0
  # pv_inj = 0.3

  C_high_sal=1.0
  C_low_sal=0.0
  SF_high_sal=0.0 # similar to oil wet
  SF_low_sal=1.0 # similar to water wet

  Nx = 500 # number of cells in x direction
  H = core_props.length # [m] length of the domain in y direction
  D_core = core_props.diameter
  A_core=π*D_core^2/4
  V_core=H*A_core
  m = createMesh1D(Nx, H)
  # m = createMesh2D(Nx, Ny, H, H)
  # m = createMeshCylindrical2D(Nx, Ny, W, H) # creates a 2D mesh
  ## define the physical parametrs
  # all the Corey-type relperm parametrs are defined for the oil-wet and water-wet
  # cases
  k0 = core_props.permeability # [m^2] average reservoir permeability
  phi0 = core_props.porosity # average porosity
  mu_oil_pw = fluids_pw.oil_viscosity # [Pa.s] oil viscosity
  mu_water_pw = fluids_pw.water_viscosity # [Pa.s] water viscosity
  mu_oil_fw = fluids_fw.oil_viscosity # [Pa.s] oil viscosity
  mu_water_fw = fluids_fw.water_viscosity # [Pa.s] water viscosity
  krw0_ww = rel_perms_pw.krw0
  krw0_ow = rel_perms_fw.krw0
  kro0_ww = rel_perms_pw.kro0
  kro0_ow = rel_perms_fw.kro0
  nw_ww = rel_perms_pw.nw
  nw_ow= rel_perms_fw.nw
  no_ww = rel_perms_pw.no
  no_ow= rel_perms_fw.no
  sor_ww=rel_perms_pw.sor
  sor_ow=rel_perms_fw.sor
  swc_ww=rel_perms_pw.swc
  swc_ow=rel_perms_fw.swc

  # SF0=(c_insitu-C_high_sal)/(C_low_sal-C_high_sal)
  # SF=createFaceVariable(m, SF0) # 1 is water wet, 0 is oil wet
  # krw0=krw0_ww*SF+krw0_ow*(1-SF)
  # kro0=kro0_ww*SF+kro0_ow*(1-SF)
  # sor=sor_ww*SF+sor_ow*(1-SF)
  # swc=swc_ww*SF+swc_ow*(1-SF)
  # no= no_ww*SF+no_ow*(1-SF)
  # nw= nw_ww*SF+nw_ow*(1-SF)

  # sor0=a*cos(teta0)^2+b*cos(teta0)+c

  p0 = core_flood.back_pressure # [bar] pressure
  u_in= core_flood.injection_velocity # [m/s] equal to 1 m/day
  sw0 = core_flood.initial_water_saturation
  # sw0(10:end-10, 10:end-10)=swc+0.2
  # sw0 = swc+0.1 # initial water saturation
  sw_in = core_flood.injected_water_saturation

  eps1=1e-5
  perm_val= k0 #field2d(Nx,Ny,k0,V_dp,clx,cly)

  k=createCellVariable(m, perm_val)
  phi=createCellVariable(m, phi0)
  
  mu_water_pw_cell= createCellVariable(m, mu_water_pw)
  mu_water_fw_cell= createCellVariable(m, mu_water_fw)
  mu_water_pw_face= arithmeticMean(mu_water_pw_cell)
  mu_water_fw_face= arithmeticMean(mu_water_fw_cell)
  k_face = geometricMean(k)
  lw = k_face/mu_water_fw_face
  lo = k_face/mu_oil_fw

  ## Define the boundaries: all fixed Sw=1, fixed pressure everywhere(?)
  BCp = createBC(m) # Neumann BC for pressure
  BCs = createBC(m) # Neumann BC for saturation
  BCc = createBC(m) # Neumann BC for salinity
  # left boundary pressure gradient
  # BCp.left.a[:]=(krw(sw_in)*lw.xvalue(1,:)+kro(sw_in)*lo.xvalue(1,:)) BCp.left.b[:]=0 BCp.left.c[:]=-u_in
  # change the right boandary to constant pressure (Dirichlet)
  BCp.left.a[:].=-k0/mu_water_pw
  BCp.left.b[:].=0
  BCp.left.c[:].=u_in
  BCp.right.a[:].=0.0
  BCp.right.b[:].=1.0
  BCp.right.c[:].=p0
  # BCp.top.a[:]=0.0
  # BCp.top.b[:]=1.0
  # BCp.top.c[:]=p0
  # BCp.bottom.a[:]=0.0
  # BCp.bottom.b[:]=1.0
  # BCp.bottom.c[:]=p0
  # change the left boundary to constant saturation (Dirichlet)
  # BCs.left.a[:]=0 BCs.left.b[:]=1 BCs.left.c[:]=1.0-sor
  # calculate the water saturation at the boundaries
  # sw_b=fzero(f1, 0.5)
  # println(sw_b)
  # for test, no capillary end effect
  # BCs.right.a[:]=0.0
  # BCs.right.b[:]=1.0
  # BCs.right.c[:]=sw_b
  BCs.left.a[:].=0.0
  BCs.left.b[:].=1.0
  BCs.left.c[:].=sw_in
  # BCs.top.a[:]=0.0
  # BCs.top.b[:]=1.0
  # BCs.top.c[:]=sw_b
  # BCs.bottom.a[:]=0.0
  # BCs.bottom.b[:]=1.0
  # BCs.bottom.c[:]=sw_b
  # Boundary conditions for salinity
  BCc.left.a[:].=0.0
  BCc.left.b[:].=1.0
  BCc.left.c[:].=c_inj

  ## define the time step and solver properties
  # dt = 1000 # [s] time step
  # dt=(W/Nx)/u_in/20 # [s]
  # pv_inj=5
  t_end = pv_inj*V_core*phi0/(u_in*A_core) # [s] final time
  dt=H/Nx/u_in/20
  SW0 = sw0
  eps_p = 1e-3 # pressure accuracy
  eps_sw = 1e-6 # saturation accuracy
  ## define the variables
  sw_old = createCellVariable(m, SW0, BCs)
  sw_new=createCellVariable(m, SW0, BCs)
  p_old = createCellVariable(m, p0, BCp)
  p_new= createCellVariable(m, p0, BCp)
  c_old = createCellVariable(m, c_insitu, BCc)
  c_face=arithmeticMean(c_old) # just to initialize
  SF=(1/(C_low_sal-C_high_sal))*(c_face-C_high_sal)
  # calculate the new relperm curves
  krw0=krw0_ww*SF+krw0_ow*(1-SF)
  kro0=kro0_ww*SF+kro0_ow*(1-SF)
  sor=sor_ww*SF+sor_ow*(1-SF)
  swc=swc_ww*SF+swc_ow*(1-SF)
  no= no_ww*SF+no_ow*(1-SF)
  nw= nw_ww*SF+nw_ow*(1-SF)

  # rel_perm_av = oil_water_rel_perms(krw0 = krw0, kro0=kro0, swc=swc, sor=sor, nw=nw, no=no)


  sw = copyCell(sw_old)
  oil_init=domainInt(1-sw_old)
  p = copyCell(p_old)
  uw = -gradientTerm(p_old) # an estimation of the water velocity
  ## start the main loop
  # generate intial pressure profile (necessary to initialize the fully
  # implicit solver)
  rec_fact=zeros(1)
  c_out_sal=zeros(1)
  water_cut=zeros(1)
  c_out_sal[1]=c_face.xvalue[end]
  t_sec=zeros(1)
  t = 0.0
  dt0=dt
  dsw_alwd= 0.005
  dp_alwd= 100.0 # Pa
  dc_alwd= 0.005
  labdaw=createFaceVariable(m,1.0) # initialize
  labdao=createFaceVariable(m,1.0) # initialize
  FL=fluxLimiter("SUPERBEE")

  prog_1=Progress(100, 1)

  while (t<t_end)
  # print("progress is $((t/t_end*100)) %", "\u1b[1G")
  ProgressMeter.update!(prog_1, floor(Int, t/t_end*100))
  # Implicit loop
  error_p = 1e5
  error_sw = 1e5
  # implicit solver
  loop_count=0
  while ((error_sw>eps_sw))
    loop_count+=1
    if loop_count>10
      break
    end
    # calculate parameters
    pgrad = gradientTerm(p)
    uw = -labdaw*pgrad
  #         pcgrad=gradientTerm(pc(sw))
    sw_face = upwindMean(sw, uw) # average value of water saturation
    # sw_ave=arithmeticMean(sw)
    # solve for pressure at known Sw
    labdao = lo*faceEval(CF.kro, sw_face, kro0, sor, swc, no)
    labdaw = lw*faceEval(CF.krw, sw_face, krw0, sor, swc, nw)
    dlabdaodsw = lo*faceEval(CF.dkrodsw, sw_face, kro0, sor, swc, no)
    dlabdawdsw = lw*faceEval(CF.dkrwdsw, sw_face, krw0, sor, swc, nw)
    labda = labdao+labdaw
    dlabdadsw = dlabdaodsw+dlabdawdsw
    # compute [Jacobian] matrices
    Mdiffp1 = diffusionTerm(-labda)
    Mdiffp2 = diffusionTerm(-labdaw)
    # Mconvsw1 = convectionUpwindTerm(-dlabdadsw.*pgrad, uw)
    # Mconvsw2 = convectionUpwindTerm(-dlabdawdsw.*pgrad, uw)
    Mconvsw1 = convectionUpwindTerm(-dlabdadsw*pgrad, uw)
    Mconvsw2 = convectionUpwindTerm(-dlabdawdsw*pgrad, uw)
    Mtranssw2, RHStrans2 = transientTerm(sw_old, dt, phi)
    # Compute RHS values
    RHS1 = Mconvsw1*sw.value # divergenceTerm(-dlabdadsw.*sw_face.*pgrad)
    RHS2 = Mconvsw2*sw.value # divergenceTerm(-dlabdawdsw.*sw_face.*pgrad)
    # include boundary conditions
    Mbcp, RHSbcp = boundaryConditionTerm(BCp)
    Mbcsw, RHSbcsw = boundaryConditionTerm(BCs)
    # Couple the equations; BC goes into the block on the main diagonal
    M = [Mdiffp1+Mbcp Mconvsw1; Mdiffp2 Mconvsw2+Mtranssw2+Mbcsw]
    RHS = [RHS1+RHSbcp; RHS2+RHStrans2+RHSbcsw]
    x = M\RHS
    p_new.value[:] .= (x[1:(Nx+2)])
    sw_new.value[:] .= (x[(Nx+2)+1:end])

    error_p = maximum(abs.((internalCells(p_new).-internalCells(p))./internalCells(p_new)))
    error_sw = maximum(abs.(internalCells(sw_new).-internalCells(sw)))
    # dt_new=dt*min(dp_alwd/error_p, dsw_alwd/error_sw)
    # print("sw_error = $error_sw \n")
    # assign new values of p and sw
    p.value[:]=p_new.value[:]
    sw.value[:]=sw_new.value[:]
  end # end of implicit loop
  if loop_count>10
    p=copyCell(p_old)
    sw=copyCell(sw_old)
    dt=dt/5
    continue
  end
  # solve the advection equation
  # ϕ d/dt(c sw)=[sw dc/dt + c dsw/dt]
  dswdt=(sw_new-sw_old)/dt
  Msc=linearSourceTerm(phi*dswdt)
  sw_face = upwindMean(sw, -labdaw*pgrad) # average value of water saturation
  pgrad = gradientTerm(p_new)
  labdaw = lw*faceEval(CF.krw, sw_face, krw0, sor, swc, nw)
  uw=-labdaw*pgrad
  Mtc, RHStc=transientTerm(c_old, dt, phi*sw)
  Mbcc, RHSbcc=boundaryConditionTerm(BCc)
  c_temp=copyCell(c_old)
  for j in 1:4
    Mconvc, RHSconv=convectionTvdTerm(uw, c_temp, FL)
    c_temp=solvePDE(m, Msc+Mtc+Mconvc+Mbcc, RHStc+RHSbcc+RHSconv)
  end
  c_new=copyCell(c_temp)
  c_face=tvdMean(c_new, uw, FL)

  # calculate the wettability modifier
  SF=(1/(C_low_sal-C_high_sal))*(c_face-C_high_sal)
  # SF.xvalue[SF.xvalue.>1.0]=1.0
  # SF.xvalue[SF.xvalue.<0.0]=0.0

  # calculate the new relperm curves
  krw0=krw0_ww*SF+krw0_ow*(1-SF)
  kro0=kro0_ww*SF+kro0_ow*(1-SF)
  sor=sor_ww*SF+sor_ow*(1-SF)
  swc=swc_ww*SF+swc_ow*(1-SF)
  no= no_ww*SF+no_ow*(1-SF)
  nw= nw_ww*SF+nw_ow*(1-SF)

  # Calculate the viscosity
  mu_w_face = SF*mu_water_pw_face+(1-SF)*mu_water_fw_face
  lw = k_face/mu_w_face

  # adaptive time step
  dp = maximum(abs.((internalCells(p_new)-internalCells(p_old))./internalCells(p_new)))
  dsw = maximum(abs.(internalCells(sw_new)-internalCells(sw_old)))
  dc = maximum(abs.(internalCells(c_new)-internalCells(c_old)))
  # dt=min(dt*(dsw_alwd/dsw), dt*(dc_alwd/dc), 15*dt, t_end-t-eps())

  # update the variables
  t+=dt
  c_old=copyCell(c_new)
  p_old = copyCell(p)
  sw_old = copyCell(sw)

  # calculate the recovery factor
  push!(rec_fact, (oil_init-domainInt(1-sw))/oil_init)
  push!(t_sec, t)
  push!(c_out_sal, c_face.xvalue[end])
  push!(water_cut, labdaw.xvalue[end]/(labdaw.xvalue[end]+labdao.xvalue[end]))

  # print(t)
  #GR.plot(sw.value[2:end-1])
  # GR.plot(c_new.value[2:end-1])
  #GR.plot(SF.xvalue)
  # GR.imshow(sw.value[2:end-1,2:end-1])
  # visualizeCells(1-sw)
  # GR.plot(t_sec/3600/24, rec_fact)
  # xlabel('time [day]')
  # ylabel('recovery factor')
  # title([num2str(t/3600/24) ' day']) drawnow
  end
  pv = t_sec*u_in/(H*phi0)
  sw_face = linearMean(sw)
  c_face  = linearMean(c_old)
  x = m.facecenters.x
  t_sec, pv, rec_fact, x, sw_face.xvalue, c_face.xvalue, c_out_sal, water_cut
end # imb_impes
