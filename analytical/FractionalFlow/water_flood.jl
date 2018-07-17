function water_flood(core_props, fluids, rel_perms, core_flood)
    fw, dfw = fractional_flow_function(rel_perms, fluids)
    sw_init = core_flood.initial_water_saturation
    sw_shock = tangent_line_saturation(rel_perms, fluids, (sw_init, fw(sw_init)))
    # println("sw_shock = $sw_shock")
    t_D_BT = 1/dfw(sw_shock) # breakthrough (BT) time [#PV]
    # println("breakthrough time = $t_D_BT")

    # tracer
    sw_tracer = tangent_line_saturation(rel_perms, fluids, (0.0, 0.0))
    t_D_tracer = 1/dfw(sw_tracer)
    xt_tracer_shock = dfw(sw_tracer)
    # construct the recovery factor curve versus the # of PV
    # R = zeros(1)
    # pv_R = zeros(1)
    # at breakthrough
    # push!(R, (1-fw(sw_init))*t_D_BT/(1-sw_init)) # recovery at BT
    # push!(pv_R, t_D_BT) # BT time
    # after breakthrough
    pv_inj = max(core_flood.injected_pore_volume, 2.0) # at least inject 2 pv
    # sw_max = outlet_saturation(pv_inj, sw_shock, dfw, rel_perms)
    # println(sw_max)
    # sw_tmp = linspace(sw_shock, sw_max, 100)
    # t_D_tmp = 1./dfw.(sw_tmp)

    # s_av_tmp = sw_tmp-(fw.(sw_tmp)-1).*t_D_tmp
    # R_tmp = (s_av_tmp-sw_init)/(1-sw_init)
    # append!(R, R_tmp[2:end])
    # append!(pv_R, t_D_tmp[2:end])

    # saturation profile
    sor = rel_perms.sor
    sw_inj = core_flood.injected_water_saturation
    phi = core_props.porosity
    ut = core_flood.injection_velocity
    L = core_props.length
    pv_to_t = phi*L/ut
    s1 = collect(linspace(min(sw_inj, 1-sor-eps()), sw_shock, 100))
    xt_s1 = dfw.(s1)
    xt_shock = dfw(sw_shock)
    xt_prf=[xt_s1; xt_shock+eps(); 2*xt_shock]
    sw_prf=[s1; sw_init; sw_init]

    xt_tracer = [0.0, xt_tracer_shock, xt_tracer_shock+eps(), xt_prf[end]]
    c_tracer = [1.0, 1.0, 0.0, 0.0]

    # calculate the pressure drop
    k = core_props.permeability
    krw, kro, dkrw, dkro = rel_perm_functions(rel_perms)
    muo, muw = fluids.oil_viscosity, fluids.water_viscosity
    
    xt = copy(xt_prf)
    sw = copy(sw_prf)
    i=1
    while(true)
        if i==length(xt)
            break
        elseif xt[i]>=xt[i+1]
            deleteat!(xt, i+1)
            deleteat!(sw, i+1)
        else
            i+=1
        end
    end
    
    
    x = collect(linspace(0,L,200))
    # println(xt)
    sw_int = Spline1D(ut/phi.*xt, sw, k=1, bc="nearest")
    t_inj=pv_inj*pv_to_t
    t = collect(linspace(0.0,t_inj, 200)) # [s] time
    p_inj = zeros(length(t))
    R_int = zeros(length(t))
    p_inj[1]=L*ut./(k*(kro.(sw_init)/muo+krw.(sw_init)/muw))
    for i in 2:length(t)
        xt_real = x/t[i]
        p_inj[i] = trapz(x, ut./(k*(kro.(sw_int(xt_real))/muo+krw.(sw_int(xt_real))/muw)))
        R_int[i] = (trapz(x, sw_int(xt_real))/L-sw_init)/(1-sw_init)
    end

    return FracFlowResults([fw], [Line([sw_init, fw(sw_init)], [sw_shock, fw(sw_shock)])],
                            [t/pv_to_t R_int], [t R_int], [xt_prf sw_prf], [xt_tracer c_tracer],
                            [t/pv_to_t p_inj], [t p_inj])
    # return pv_R, R, xt_prf, sw_prf # for the time being to test the code

end
