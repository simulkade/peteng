

function single_ion_adsorption_water_flood_single_shock(core_props, fluids_ls, fluids_hs, rel_perms_hs,
    rel_perms_ls, core_flood, eq_const)
    # ====================================================================
    # This function is not done yet. The trick is to draw a tangent from
    # the high sal 1-sor_hs to the low sal curve to find the shock front
    # saturation. The low sal shock speed is calculated at this saturation
    # by fw/(sw+eq_const)
    # The implementation needs some if statements to make sure that the
    # shock front saturation calculated by the tangent line from -eq_const
    # to the low sal fw is indeed lower than 1-sor_hs
    # ====================================================================
    # construct the fractional flow curves
    fw_ls, dfw_ls = fractional_flow_function(rel_perms_ls, fluids_ls)
    fw_hs, dfw_hs = fractional_flow_function(rel_perms_hs, fluids_hs)
    # low sal shock (tangent line from (0,0))
    sw_shock_ls = tangent_line_saturation(rel_perms_ls, fluids_ls, (-eq_const, 0.0))
    # println("low sal sw_shock = $sw_shock_ls")
    t_D_BT_ls = 1/dfw_ls(sw_shock_ls) # breakthrough (BT) time [#PV]
    # println("low sal breakthrough time = $t_D_BT_ls")
    # High sal shock (cross point between the ls tangent and the hs fw)
    sw_shock_hs = cross_point_saturation(fw_hs, rel_perms_hs, (-eq_const, 0.0),
            (sw_shock_ls, fw_ls(sw_shock_ls)))
    # println("high sal sw_shock = $sw_shock_hs")
    sw_init = core_flood.initial_water_saturation
    t_D_BT_hs = (sw_shock_hs-sw_init)/(fw_hs(sw_shock_hs)-fw_hs(sw_init)) # breakthrough (BT) time [#PV]
    # println("high sal breakthrough time = $t_D_BT_hs")
    if t_D_BT_hs<t_D_BT_ls
        info("Using the single_ion_adsorption_water_flood function.")
        return single_ion_adsorption_water_flood(core_props, fluids_ls, fluids_hs, rel_perms_hs,
                rel_perms_ls, core_flood, eq_const)
    end
    sor_ls = rel_perms_ls.sor
    sw_shock = cross_point_saturation(fw_ls, rel_perms_ls, (sw_shock_hs, fw_hs(sw_shock_hs)),
        (sw_init, fw_hs(sw_init)), sw_left = sw_shock_ls, sw_right = 1-sor_ls)
    # sw_shock = cross_point_saturation(fw_ls, rel_perms_ls, (-eq_const, 0.0),
    #         (sw_init, fw_hs(sw_init)), sw_left = sw_shock_ls, sw_right = 1-sor_ls) # just a test
    # sw_shock = cross_point_saturation(fw_ls, rel_perms_ls, (sw_shock_ls, fw_ls(sw_shock_ls)),
    #     (sw_init, fw_hs(sw_init)), sw_left = sw_shock_ls, sw_right = 1-sor_ls)
    if sw_shock<sw_shock_ls
        sw_shock = cross_point_saturation(fw_ls, rel_perms_ls, (sw_shock_hs, fw_hs(sw_shock_hs)),
            (sw_init, fw_hs(sw_init)), sw_left = sw_shock_ls, sw_right = 1-sor_ls)
    end
    sor_hs = rel_perms_hs.sor
    # t_D_BT_hs = (sw_shock-sw_init)/(fw_ls(sw_shock)-fw_hs(sw_init))
    t_D_BT_hs = (sw_shock+eq_const)/(fw_ls(sw_shock)-0.0)
    # println("sw_shock = $sw_shock")
    t_D_BT_ls = 1/dfw_ls(sw_shock)
    # println("low sal breakthrough time = $t_D_BT_ls")

    # tracer
    # I'm not sure about this, but the numerical solution suggests that
    # the tracer front follows the shock front for this special case
    # sw_tracer = tangent_line_saturation(rel_perms_ls, fluids_ls, (0.0, 0.0))
    # sw_tracer = max(sw_tracer, sw_shock)
    t_D_tracer = t_D_BT_hs # (sw_shock-sw_init)/(fw_ls(sw_shock)-fw_hs(sw_init))
    xt_tracer_shock = 1/t_D_tracer

    # construct the recovery factor curve versus the # of PV
    # R = zeros(1)
    # pv_R = zeros(1)
    # at breakthrough of the hs brine
    sw_init = core_flood.initial_water_saturation
    # push!(R, 0.0) # recovery at BT
    # push!(pv_R, t_D_BT_ls) # BT time
    # at breakthrough of the ls brine
    # push!(R, (1-fw_hs(sw_shock))*(t_D_BT_ls-t_D_BT_hs)/(1-sw_init)) # recovery at BT
    # push!(pv_R, t_D_BT_ls) # BT time

    # after breakthrough
    pv_inj = max(core_flood.injected_pore_volume, 2.0) # at least inject 2 pv
    # f_sw = sw -> (pv_inj-1/dfw_ls(sw)) # find the outlet saturation at the end of injection
    # sw_max = fzero(f_sw, sw_shock)
    # println("swmax $sw_max")
    # sw_tmp = linspace(sw_shock, sw_max, 100)
    # t_D_tmp = 1./dfw_ls.(sw_tmp)

    # s_av_tmp = sw_tmp-(fw_ls.(sw_tmp)-1).*t_D_tmp
    # R_tmp = (s_av_tmp-sw_init)/(1-sw_init)
    # append!(R, R_tmp)
    # append!(pv_R, t_D_tmp)

    # saturation profile
    sw_inj = core_flood.injected_water_saturation
    phi = core_props.porosity
    ut = core_flood.injection_velocity
    L = core_props.length
    pv_to_t = phi*L/ut
    s1 = collect(linspace(min(sw_inj, 1-sor_ls-eps()), sw_shock, 100))
    xt_s1 = dfw_ls.(s1)
    xt_shock_ls = 1/t_D_BT_ls
    xt_shock_hs = 1/t_D_BT_hs
    xt_prf=[xt_s1; xt_shock_ls+eps(); xt_shock_hs; xt_shock_hs+eps(); 1/0.3]
    sw_prf=[s1; sw_shock; sw_shock; sw_init; sw_init]
    xt_tracer = [0.0, xt_tracer_shock, xt_tracer_shock+eps(), xt_prf[end]]
    c_tracer = [1.0, 1.0, 0.0, 0.0]

    # extract data for calculating the pressure drop
    k = core_props.permeability
    krw_ls, kro_ls, dkrw_ls, dkro_ls = rel_perm_functions(rel_perms_ls)
    muo_ls, muw_ls = fluids_ls.oil_viscosity, fluids_ls.water_viscosity
    krw_hs, kro_hs, dkrw_hs, dkro_hs = rel_perm_functions(rel_perms_hs)
    muo_hs, muw_hs = fluids_hs.oil_viscosity, fluids_hs.water_viscosity
    λ_ls = s -> krw_ls(s)/muw_ls+kro_ls(s)/muo_ls
    λ_hs = s -> krw_hs(s)/muw_hs+kro_hs(s)/muo_hs
    λ_t = [λ_ls.(s1); λ_hs(sw_shock); λ_hs(sw_shock); λ_hs(sw_init); λ_hs(sw_init)]

    xt = copy(xt_prf)
    sw = copy(sw_prf)
    i=1
    while(true)
        if i==length(xt)
            break
        elseif xt[i]>=xt[i+1]
            deleteat!(xt, i+1)
            deleteat!(sw, i+1)
            deleteat!(λ_t, i+1)
        else
            i+=1
        end
    end

    x = collect(linspace(0,L,200))
    sw_int = Spline1D(ut/phi.*xt, sw, k=1, bc="nearest")
    λ_int = Spline1D(ut/phi.*xt, λ_t, k=1, bc="nearest")
    t_inj=pv_inj*pv_to_t
    t = collect([linspace(0,t_D_BT_hs*pv_to_t, 10); 
                linspace((t_D_BT_hs+eps())*pv_to_t, t_D_BT_ls*pv_to_t, 40);
                linspace((t_D_BT_ls+eps())*pv_to_t,t_inj, 100)]) # [s] time
    p_inj = zeros(length(t))
    R_int = zeros(length(t))
    p_inj[1]=L*ut./(k*λ_hs(sw_init))
    for i in 2:length(t)
        xt_real = x/t[i]
        p_inj[i] = trapz(x, ut./(k*λ_int(xt_real)))
        R_int[i] = (trapz(x, sw_int(xt_real))/L-sw_init)/(1-sw_init)
    end

    R_int[R_int.<0.0] = 0.0
 

    return FracFlowResults([fw_hs, fw_ls], [Line([-eq_const, 0.0], [sw_shock_ls, fw_ls(sw_shock_ls)]),
                            Line([sw_init, fw_hs(sw_init)], [sw_shock_hs, fw_hs(sw_shock_hs)]),
                            Line([sw_shock, fw_ls(sw_shock)], [sw_init, fw_hs(sw_init)])],
                            [t/pv_to_t R_int], [t R_int], [xt sw], [xt_tracer c_tracer],
                            [t/pv_to_t p_inj], [t p_inj])
    # return pv_R, R, xt_prf, sw_prf # for the time being to test the code
end

