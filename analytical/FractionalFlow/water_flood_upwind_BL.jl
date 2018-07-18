"""
Water flooding with BL assumptions
"""
function water_flood_numeric(core_props, fluids, rel_perms, core_flood; Nx = 30)
    pv_inj = core_flood.injected_pore_volume

    # Nx = 100 # number of cells in x direction
    H = core_props.length # [m] length of the domain in y direction
    D_core = core_props.diameter
    A_core=π*D_core^2/4
    V_core=H*A_core
    m = createMesh1D(Nx, H)
    
    ## define the physical parametrs
    # all the Corey-type relperm parametrs are defined for the oil-wet and water-wet
    # cases
    k0 = core_props.permeability # [m^2] average reservoir permeability
    phi0 = core_props.porosity # average porosity
    mu_oil = fluids.oil_viscosity # [Pa.s] oil viscosity
    mu_water = fluids.water_viscosity # [Pa.s] water viscosity
    krw0 = rel_perms.krw0
    kro0 = rel_perms.kro0
    nw = rel_perms.nw
    no = rel_perms.no
    sor=rel_perms.sor
    swc=rel_perms.swc

    p0 = core_flood.back_pressure # [bar] pressure
    u_in= core_flood.injection_velocity # [m/s] equal to 1 m/day
    sw0 = core_flood.initial_water_saturation

    sw_in = core_flood.injected_water_saturation

    eps1=1e-5
    perm_val= k0 #field2d(Nx,Ny,k0,V_dp,clx,cly)

    k=createCellVariable(m, perm_val)
    phi=createCellVariable(m, phi0)

    lw = geometricMean(k)/mu_water
    lo = geometricMean(k)/mu_oil

    ## Define the boundaries: all fixed Sw=1, fixed pressure everywhere(?)
    BCp = createBC(m) # Neumann BC for pressure
    BCs = createBC(m) # Neumann BC for saturation
    BCc = createBC(m) # Neumann BC for salinity

    BCp.left.a[:]=-k0/mu_water
    BCp.left.b[:]=0
    BCp.left.c[:]=u_in

    BCp.right.a[:]=0.0
    BCp.right.b[:]=1.0
    BCp.right.c[:]=p0

    BCs.left.a[:]=0.0
    BCs.left.b[:]=1.0
    BCs.left.c[:]=sw_in

    BCc.left.a[:]=0.0
    BCc.left.b[:]=1.0
    BCc.left.c[:]=1.0 # injection concentration

    ## define the time step and solver properties
    t_end = pv_inj*V_core*phi0/(u_in*A_core) # [s] final time
    dt=H/Nx/u_in/10
    SW0 = sw0
    eps_p = 1e-3 # pressure accuracy
    eps_sw = 1e-6 # saturation accuracy
    ## define the variables
    sw_old = createCellVariable(m, SW0, BCs)
    sw_new=createCellVariable(m, SW0, BCs)
    p_old = createCellVariable(m, p0, BCp)
    p_new= createCellVariable(m, p0, BCp)
    c_old = createCellVariable(m, 0.0, BCc)
    c_face=arithmeticMean(c_old) # just to initialize

    # define the rel perm functions
    KRW = sw -> CF.krw(sw, krw0, sor, swc, nw)
    KRO = sw -> CF.kro(sw, kro0, sor, swc, no)
    dKRW = sw -> CF.dkrwdsw(sw, krw0, sor, swc, nw)
    dKRO = sw -> CF.dkrodsw(sw, kro0, sor, swc, no)

    sw = copyCell(sw_old)
    oil_init=domainInt(1-sw_old)
    p = copyCell(p_old)
    uw = -gradientTerm(p_old) # an estimation of the water velocity
    ## start the main loop
    # generate intial pressure profile (necessary to initialize the fully
    # implicit solver)
    dp_core = zeros(1)
    dp_core[1]=H*u_in/(k0*(CF.kro(sw0, kro0, sor, swc, no)/mu_oil+
        CF.krw(sw0, krw0, sor, swc, nw)/mu_water))
    rec_fact=zeros(1)
    c_out_sal=zeros(1)
    c_out_sal[1]=c_face.xvalue[end]
    t_sec=zeros(1)
    t = 0.0
    dt0=dt
    dsw_alwd= 0.005
    dp_alwd= 100.0 # Pa
    dc_alwd= 0.005
    labdaw=createFaceVariable(m,1.0) # initialize
    FL=fluxLimiter("SUPERBEE")

    prog_1=Progress(100, 1)

    while (t<t_end)
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
            uw = -labdaw.*pgrad
            sw_face = upwindMean(sw, uw) # average value of water saturation
            # sw_ave=arithmeticMean(sw)
            # solve for pressure at known Sw
            labdao = lo.*faceEval(KRO, sw_face)
            labdaw = lw.*faceEval(KRW, sw_face)
            dlabdaodsw = lo.*faceEval(dKRO, sw_face)
            dlabdawdsw = lw.*faceEval(dKRW, sw_face)
            labda = labdao+labdaw
            dlabdadsw = dlabdaodsw+dlabdawdsw
            # compute [Jacobian] matrices
            Mdiffp1 = diffusionTerm(-labda)
            Mdiffp2 = diffusionTerm(-labdaw)
            # Mconvsw1 = convectionUpwindTerm(-dlabdadsw.*pgrad, uw)
            # Mconvsw2 = convectionUpwindTerm(-dlabdawdsw.*pgrad, uw)
            Mconvsw1 = convectionUpwindTerm(-dlabdadsw.*pgrad, uw)
            Mconvsw2 = convectionUpwindTerm(-dlabdawdsw.*pgrad, uw)
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
            p_new.value[:] = full(x[1:(Nx+2)])
            sw_new.value[:] = full(x[(Nx+2)+1:end])

            error_p = maximum(abs.((internalCells(p_new)-internalCells(p))./internalCells(p_new)))
            error_sw = maximum(abs.(internalCells(sw_new)-internalCells(sw)))
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
        Msc=linearSourceTerm(phi.*dswdt)
        sw_face = upwindMean(sw, -labdaw.*pgrad) # average value of water saturation
        pgrad = gradientTerm(p_new)
        labdaw = lw.*faceEval(KRW, sw_face)
        uw=-labdaw.*pgrad
        Mtc, RHStc=transientTerm(c_old, dt, phi.*sw)
        Mbcc, RHSbcc=boundaryConditionTerm(BCc)
        c_temp=copyCell(c_old)
        for j in 1:4
            Mconvc, RHSconv=convectionTvdTerm(uw, c_temp, FL)
            c_temp=solveLinearPDE(m, Msc+Mtc+Mconvc+Mbcc, RHStc+RHSbcc+RHSconv)
        end
        c_new=copyCell(c_temp)
        c_face=tvdMean(c_new, uw, FL)

        # adaptive time step
        dp = maximum(abs.((internalCells(p_new)-internalCells(p_old))./internalCells(p_new)))
        dsw = maximum(abs.(internalCells(sw_new)-internalCells(sw_old)))
        dc = maximum(abs.(internalCells(c_new)-internalCells(c_old)))

        # update the variables
        t+=dt
        update!(prog_1, floor(Int, t/t_end*100))
        c_old=copyCell(c_new)
        p_old = copyCell(p)
        sw_old = copyCell(sw)

        dt=min(dt*(dsw_alwd/dsw), 15*dt, t_end-t-eps())

        # calculate the recovery factor
        push!(rec_fact, (oil_init-domainInt(1-sw))/oil_init)
        push!(t_sec, t)
        push!(c_out_sal, c_face.xvalue[end])
        push!(dp_core, 0.5*sum(p.value[1:2])-p0)

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
    t_sec, pv, rec_fact, dp_core, x, sw_face.xvalue, c_face.xvalue, c_out_sal
end # imb_impes
