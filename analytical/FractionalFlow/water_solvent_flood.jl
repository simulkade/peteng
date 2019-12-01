"""
the eq_const is defined as the ratio of the volume fraction of the solvent in oil to
its volume fraction in water. For instance, For DME, the molar fraction is around 2.0
One criteria is that eq_const should probably be higher than 1.
"""
function water_soluble_solvent_flood(core_props, fluids_ls, fluids_hs, rel_perms_hs,
    rel_perms_ls, core_flood, eq_const)
    # construct the fractional flow curves
    fw_ls, dfw_ls = fractional_flow_function(rel_perms_ls, fluids_ls)
    fw_hs, dfw_hs = fractional_flow_function(rel_perms_hs, fluids_hs)
    # low sal shock (tangent line from (0,0))
    point1 = (-eq_const/(1-eq_const), -eq_const/(1-eq_const))
    sw_shock_ls = tangent_line_saturation(rel_perms_ls, fluids_ls, point1)
    println("low sal sw_shock = $sw_shock_ls")
    t_D_BT_ls = 1/dfw_ls(sw_shock_ls) # breakthrough (BT) time [#PV]
    println("low sal breakthrough time = $t_D_BT_ls")
    # High sal shock (cross point between the ls tangent and the hs fw)
    sw_shock_hs = cross_point_saturation(fw_hs, rel_perms_hs, point1, (sw_shock_ls, fw_ls(sw_shock_ls)))
    # println("high sal sw_shock = $sw_shock_hs")
    sw_init = core_flood.initial_water_saturation
    t_D_BT_hs = (sw_shock_hs-sw_init)/(fw_hs(sw_shock_hs)-fw_hs(sw_init)) # breakthrough (BT) time [#PV]
    # println("high sal breakthrough time = $t_D_BT_hs")

    # tracer
    sw_tracer = tangent_line_saturation(rel_perms_ls, fluids_ls, (0.0, 0.0))
    t_D_tracer = 1/dfw_ls(sw_tracer)
    xt_tracer_shock = dfw_ls(sw_tracer)

    # construct the recovery factor curve versus the # of PV
    # R = zeros(1)
    # pv_R = zeros(1)
    # at breakthrough of the hs brine
    # push!(R, (1-fw_hs(sw_init))*t_D_BT_hs/(1-sw_init)) # recovery at BT
    # push!(pv_R, t_D_BT_hs) # BT time
    # at breakthrough of the ls brine
    # push!(R, R[end]+(1-fw_hs(sw_shock_hs))*(t_D_BT_ls-t_D_BT_hs)/(1-sw_init)) # recovery at BT
    # push!(pv_R, t_D_BT_ls) # BT time

    # after breakthrough
    pv_inj = max(core_flood.injected_pore_volume, 2.0) # at least inject 2 pv
    # f_sw = sw -> (pv_inj-1/dfw_ls(sw)) # find the outlet saturation at the end of injection
    # sw_max = fzero(f_sw, sw_shock_ls)
    # sw_tmp = linspace(sw_shock_ls, sw_max, 100)
    # t_D_tmp = 1.0./dfw_ls.(sw_tmp)

    # s_av_tmp = sw_tmp-(fw_ls.(sw_tmp)-1).*t_D_tmp
    # R_tmp = (s_av_tmp-sw_init)/(1-sw_init)
    # append!(R, R_tmp)
    # append!(pv_R, t_D_tmp)

    # saturation profile
    sor = rel_perms_ls.sor
    sw_inj = core_flood.injected_water_saturation
    phi = core_props.porosity
    ut = core_flood.injection_velocity
    L = core_props.length
    pv_to_t = phi*L/ut
    s1 = collect(linspace(min(sw_inj, 1-sor-eps()), sw_shock_ls, 300))
    xt_s1 = dfw_ls.(s1)
    xt_shock_ls = 1/t_D_BT_ls
    xt_shock_hs = 1/t_D_BT_hs
    xt_prf=[xt_s1; xt_shock_ls+eps(); xt_shock_hs; xt_shock_hs+eps(); 1/0.3]
    sw_prf=[s1; sw_shock_hs; sw_shock_hs; sw_init; sw_init]
    xt_tracer = [0.0, xt_tracer_shock, xt_tracer_shock+eps(), xt_prf[end]]
    c_tracer = [1.0, 1.0, 0.0, 0.0]
    # return pv_R, R, xt_prf, sw_prf # for the time being to test the code

    # extract data for calculating the pressure drop
    k = core_props.permeability
    krw_ls, kro_ls, dkrw_ls, dkro_ls = rel_perm_functions(rel_perms_ls)
    muo_ls, muw_ls = fluids_ls.oil_viscosity, fluids_ls.water_viscosity
    krw_hs, kro_hs, dkrw_hs, dkro_hs = rel_perm_functions(rel_perms_hs)
    muo_hs, muw_hs = fluids_hs.oil_viscosity, fluids_hs.water_viscosity
    λ_ls = s -> krw_ls(s)/muw_ls+kro_ls(s)/muo_ls
    λ_hs = s -> krw_hs(s)/muw_hs+kro_hs(s)/muo_hs
    λ_w_ls = s -> krw_ls(s)/muw_ls
    λ_w_hs = s -> krw_hs(s)/muw_hs
    λ_t = [λ_ls.(s1); λ_hs(sw_shock_hs); λ_hs(sw_shock_hs); λ_hs(sw_init); λ_hs(sw_init)]
    λ_w = [λ_w_ls.(s1); λ_w_hs(sw_shock_hs); λ_w_hs(sw_shock_hs); λ_w_hs(sw_init); λ_w_hs(sw_init)]

    xt = copy(xt_prf)
    sw = copy(sw_prf)
    i=1
    while(true)
        if i==length(xt)
            break
        elseif xt[i]>=xt[i+1]
            deleteat!(xt, i)
            deleteat!(sw, i)
            deleteat!(λ_t, i)
            deleteat!(λ_w, i)
        else
            i+=1
        end
    end

    x = collect(linspace(0,L,200))
    sw_int = Spline1D(xt, sw, k=1, bc="nearest")
    λ_int = Spline1D(xt, λ_t, k=1, bc="nearest")
    λ_w_int = Spline1D(xt, λ_w, k=1, bc="nearest")
    t_inj=pv_inj*pv_to_t
    t = collect([linspace(0,t_D_BT_hs*pv_to_t, 10);
                linspace((t_D_BT_hs+eps())*pv_to_t, t_D_BT_ls*pv_to_t, 40);
                linspace((t_D_BT_ls+eps())*pv_to_t,t_inj, 100)]) # [s] time
    println("co2 breakthrough time is $(t_D_BT_ls*pv_to_t)")
    p_inj = zeros(length(t))
    R_int = zeros(length(t))
    WC_int = zeros(length(t))
    p_inj[1]=L*ut./(k*λ_hs(sw_init))
    WC_int[1]=λ_w_hs(sw_init)/λ_hs(sw_init)
    for i in 2:length(t)
        xt_real = x/t[i]/(ut/phi)
        p_inj[i] = trapz(x, ut./(k*λ_int(xt_real)))
        R_int[i] = (trapz(x, sw_int(xt_real))/L-sw_init)/(1-sw_init)
        WC_int[i] = λ_w_int(xt_real[end])/λ_int(xt_real[end])
    end

    R_int[R_int.<0.0] .= 0.0

    return FracFlowResults([fw_hs, fw_ls], [Line([-eq_const/(1-eq_const), -eq_const/(1-eq_const)],
                            [sw_shock_hs, fw_hs(sw_shock_hs)]),
                            Line([sw_init, fw_hs(sw_init)], [sw_shock_hs, fw_hs(sw_shock_hs)])],
                            [t/pv_to_t R_int], [t R_int], [xt_prf sw_prf], [xt_tracer c_tracer],
                            [t/pv_to_t p_inj], [t p_inj], [t/pv_to_t WC_int], [t WC_int])

end
