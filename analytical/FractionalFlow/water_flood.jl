function water_flood(core_props, fluids, rel_perms, core_flood)
    fw, dfw = fractional_flow_function(rel_perms, fluids)
    sw_init = core_flood.initial_water_saturation
    sw_shock = tangent_line_saturation(rel_perms, fluids, (sw_init, fw(sw_init)))
    println("sw_shock = $sw_shock")
    t_D_BT = 1/dfw(sw_shock) # breakthrough (BT) time [#PV]
    println("breakthrough time = $t_D_BT")

    # tracer
    sw_tracer = tangent_line_saturation(rel_perms, fluids, (0.0, 0.0))
    t_D_tracer = 1/dfw(sw_tracer)
    xt_tracer_shock = dfw(sw_tracer)
    # construct the recovery factor curve versus the # of PV
    R = zeros(1)
    pv_R = zeros(1)
    # at breakthrough
    push!(R, (1-fw(sw_init))*t_D_BT/(1-sw_init)) # recovery at BT
    push!(pv_R, t_D_BT) # BT time
    # after breakthrough
    pv_inj = max(core_flood.injected_pore_volume, 2.0) # at least inject 2 pv
    f_sw = sw -> (pv_inj-1/dfw(sw)) # find the outlet saturation at the end of injection
    sw_max = fzero(f_sw, sw_shock)
    sw_tmp = linspace(sw_shock, sw_max, 100)
    t_D_tmp = 1./dfw.(sw_tmp)

    s_av_tmp = sw_tmp-(fw.(sw_tmp)-1).*t_D_tmp
    R_tmp = (s_av_tmp-sw_init)/(1-sw_init)
    append!(R, R_tmp)
    append!(pv_R, t_D_tmp)

    # saturation profile
    sor = rel_perms.sor
    sw_inj = core_flood.injected_water_saturation
    phi = core_props.porosity
    ut = core_flood.injection_velocity
    L = core_props.length
    pv_to_t = phi*L/ut
    s1 = collect(linspace(min(sw_inj, 1-sor-eps()), sw_shock-eps(), 100))
    xt_s1 = dfw.(s1)
    xt_shock = dfw(sw_shock)
    xt_prf=[xt_s1; xt_shock; xt_shock+eps(); 2*xt_shock]
    sw_prf=[s1; sw_shock; sw_init; sw_init]

    xt_tracer = [0.0, xt_tracer_shock, xt_tracer_shock+eps(), xt_prf[end]]
    c_tracer = [1.0, 1.0, 0.0, 0.0]

    # calculate the pressure drop
    k = core_props.permeability
    krw, kro, dkrw, dkro = rel_perm_functions(rel_perms)

    sw_int = Spline1D(xt_prf, sw_prf, k=1)
    t_inj=pv_inj*pv_to_t
    t = collect(linspace(0.0,t_inj, 200)) # [s] time
    p_inj = zeros(length(t))
    R_oil= zeros(length(t))
    p_inj[1]=trapz(x, ut./(k*(kro.(sw0*ones(size(x)))/muo+krw_new.(sw0*ones(size(x)))/muw)))
    for i in 2:length(t)
        xt = x/t[i]
        p_inj[i] = trapz(x, ut./(k*(kro_new.(sw_int(xt))/muo+krw_new.(sw_int(xt))/muw)))
        R_oil[i]=1.0-trapz(x/L, 1.0-sw_int(xt))/(1-sw0)
    end

    return FracFlowResults([fw], [Line([sw_init, fw(sw_init)], [sw_shock, fw(sw_shock)])],
                            [pv_R R], [pv_R*pv_to_t R], [xt_prf sw_prf], [xt_tracer c_tracer])
    # return pv_R, R, xt_prf, sw_prf # for the time being to test the code

end
