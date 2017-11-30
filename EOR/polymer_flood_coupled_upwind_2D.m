% include('../functions/rel_perms_real.jl')
clc

a1 = 4.0;
a2 = 0.0;
a3 = 6.0;
sp = -0.24;
mu_w = 0.001; % Pa.s
polymer_viscosity = @(cs, cp)(mu_w*(1.0+cp.*cs.^sp.*(a1+a2*cp+a3*cp.*cp)));

dpolymer_viscosity_dcs = @(cs, cp)(sp*mu_w*cp.*cs.^(sp-1).*(a1+a2*cp+a3*cp.*cp));

dpolymer_viscosity_dcp = @(cs, cp)(mu_w*cs.^sp.*(a1+2*a2*cp+3*a3*cp.*cp));

% domain
Lx = 1.0;
Ly = 1.0;
Nx = 50;
Ny = 50;
m = createMesh2D(Nx, Ny, Lx, Ly);

% (petro)physical properties
mu_oil = 10e-3; % Pa.s
perm = 0.1e-12;
poros = 0.3;

% rel-perms
krw0 = 0.2;
kro0 = 0.8;
no  = 2.0;
nw  = 2.0;
swc  = 0.08;
sor = 0.2;

sws=@(sw)((sw>swc).*(sw<1-sor).*(sw-swc)/(1-sor-swc)+(sw>=1-sor).*ones(size(sw)));
KRO=@(sw)((sw>=swc).*kro0.*(1-sws(sw)).^no+(sw<swc).*(1+(kro0-1)/swc*sw));
KRW=@(sw)((sw<=1-sor).*krw0.*sws(sw).^nw+(sw>1-sor).*(-(1-krw0)/sor.*(1.0-sw)+1.0));
dKRW=@(sw)((sw<=1-sor).*nw.*krw0.*(1/(1-sor-swc)).*sws(sw).^(nw-1)+(sw>1-sor)*((1-krw0)/sor));
dKRO=@(sw)((sw>=swc).*(-kro0*no*(1-sws(sw)).^(no-1))/(-swc-sor+1)+(sw<swc).*((kro0-1)/swc));

% initial and boundary conditions
p_res = 100e5; % Pa reservoir (and right boundary pressure)
sw_0 = swc; % inital water saturation
sw_inj = 1.0; % injection water saturation
cp_0 = 0.0; % initial polymer concentration kg/m^3
cp_inj = 1.0; % injected polymer concentration
cs_0 = 10.0;  % initial salt concentration kg/m^3
cs_inj = 5.0; % injected salt concentration


phi = createCellVariable(m, poros*(ones(Nx,Ny)+rand(Nx,Ny)/100)); % porosity
k = createCellVariable(m, perm); % m^2 permeability
k_face   = harmonicMean(k); % harmonic averaging for perm
u_inj = 1e-5; % m/s injection Darcy velocity
pv_inj = 0.3; % injected pore volumes
t_final = pv_inj*Lx/(u_inj/poros); % s final time
dt0 = t_final/Nx/2; % s time step
dt  = t_final/Nx/2; % variable time step

BCp = createBC(m);   % pressure boundary
BCp.left.a(:) = perm/polymer_viscosity(cs_inj, cp_inj);
BCp.left.b(:) = 0.0;
BCp.left.c(:) = -u_inj;
BCp.right.a(:) = 0.0;
BCp.right.b(:) = 1.0;
BCp.right.c(:) = p_res;
[M_bc_p, RHS_bc_p] = boundaryCondition(BCp);

% saturation boundary
BCs = createBC(m);
BCs.left.a(:) = 0.0;
BCs.left.b(:) = 1.0;
BCs.left.c(:) = sw_inj; % injected water saturation
[M_bc_s, RHS_bc_s] = boundaryCondition(BCs);

% salt concentration boundary
BCcs = createBC(m);
BCcs.left.a(:) = 0.0;
BCcs.left.b(:) = 1.0;
BCcs.left.c(:) = cs_inj; % injected salt concentration
[M_bc_cs, RHS_bc_cs] = boundaryCondition(BCcs);

% polymer concentration boundary
BCcp = createBC(m);
BCcp.left.a(:) = 0.0;
BCcp.left.b(:) = 1.0;
BCcp.left.c(:) = cp_inj; % injected polymer concentration
[M_bc_cp, RHS_bc_cp] = boundaryCondition(BCcp);

% initial conditions
p_init  = createCellVariable(m, p_res);
sw_init = createCellVariable(m, sw_0);
cs_init = createCellVariable(m, cs_0);
cp_init = createCellVariable(m, cp_0);
p_val  = createCellVariable(m, p_res);
sw_val = createCellVariable(m, sw_0);
cs_val = createCellVariable(m, cs_0);
cp_val = createCellVariable(m, cp_0);

uw       = -gradientTerm(p_val);
% M_bc_p, RHS_bc_p = boundaryConditionTerm(BCp)
% M_bc_s, RHS_bc_s = boundaryConditionTerm(BCs)
% M_bc_c, RHS_bc_c = boundaryConditionTerm(BCcs)
% M_bc_t, RHS_bc_t = boundaryConditionTerm(BCt)
tol_s = 1e-7;
tol_c = 1e-7;
max_int_loop = 40;
t = 0.0;
FL = fluxLimiter('SUPERBEE'); % flux limiter

while t<t_final
    error_s = 2*tol_s;
    error_cp = 2*tol_c;
    loop_countor = 0;
    while (error_s>tol_s) || (error_cp>tol_c)
        loop_countor = loop_countor + 1;
        if loop_countor > max_int_loop
            sw_val = sw_init;
            cp_val  = cp_init;
            cs_val  = cs_init;
            p_val  = p_init;
            dt = dt/3.0;
            break
        end
        grad_p0      = gradientTerm(p_val);
        grad_s0      = gradientTerm(sw_val);

        cs_face   = upwindMean(cs_val, uw);
        cp_face   = upwindMean(cp_val, uw);
        sw_face       = upwindMean(sw_val, uw);

        mu_water_face = polymer_viscosity(cs_face, cp_face);
        d_mu_cs_face = dpolymer_viscosity_dcs(cs_face, cp_face);
        d_mu_cp_face = dpolymer_viscosity_dcp(cs_face, cp_face);
        d_cs_mu_dcs = (mu_water_face-d_mu_cs_face.*cs_face)./(mu_water_face.*mu_water_face);
        d_cs_mu_dcp = (-d_mu_cp_face.*cs_face)./(mu_water_face.*mu_water_face);
        d_cp_mu_dcp = (mu_water_face-d_mu_cp_face.*cp_face)./(mu_water_face.*mu_water_face);
        d_cp_mu_dcs = (-d_mu_cs_face.*cp_face)./(mu_water_face.*mu_water_face);
        d_mu_dcs = -d_mu_cs_face./(mu_water_face.*mu_water_face);
        d_mu_dcp = -d_mu_cp_face./(mu_water_face.*mu_water_face);

        % println('visc')
        krw_face      = KRW(sw_face);
        kro_face      = KRO(sw_face);

        labda_w_face      = k_face.*krw_face./mu_water_face;
        labda_o_face      = k_face.*kro_face/mu_oil;
        uw            = -labda_w_face.*grad_p0;
        uo            = -labda_o_face.*grad_p0;
        ut            = uw + uo;

        dlabda_w_face     = k_face.*dKRW(sw_face)./mu_water_face;
        dlabda_o_face     = k_face.*dKRO(sw_face)./mu_oil;

        % water mass balance (wmb)
        [M_t_s_wmb, RHS_t_s_wmb] = transientTerm(sw_init, dt, phi);
        M_d_p_wmb = diffusionTerm(-labda_w_face);
        M_a_s_wmb = convectionUpwindTerm(-dlabda_w_face.*grad_p0, ut);
        M_a_cs_wmb = convectionUpwindTerm(-k_face.*krw_face.*d_mu_dcs.*grad_p0, ut);
        M_a_cp_wmb = convectionUpwindTerm(-k_face.*krw_face.*d_mu_dcp.*grad_p0, ut);

        % oil mass balance (omb)
        [M_t_s_omb, RHS_t_s_omb] = transientTerm(sw_init, dt, -phi);
        M_d_p_omb = diffusionTerm(-labda_o_face);
        M_a_s_omb = convectionUpwindTerm(-(dlabda_o_face.*grad_p0), ut);

        % salt mass balance
        [M_t_cs_smb, RHS_t_cs_smb] = transientTerm(cs_init, dt, sw_val.*phi);
        [M_t_s_smb, RHS_t_s_smb] = transientTerm(sw_init, dt, cs_val.*phi);
        % M_a_cs_smb, RHS_a_cs_smb = convectionTvdTerm(uw, cs_val, FL, ut)
        M_a_s_smb = convectionUpwindTerm(-dlabda_w_face.*cs_face.*grad_p0, ut);        
        M_d_p_smb = diffusionTerm(-labda_w_face.*cs_face);
        M_a_cs_smb = convectionUpwindTerm(-k_face.*krw_face.*d_cs_mu_dcs.*grad_p0, ut);
        M_a_cp_smb = convectionUpwindTerm(-k_face.*krw_face.*d_cs_mu_dcp.*grad_p0, ut);

        % polymer mass balance
        [M_t_cp_pmb, RHS_t_cp_pmb] = transientTerm(cp_init, dt, sw_val.*phi);
        [M_t_s_pmb, RHS_t_s_pmb] = transientTerm(sw_init, dt, cp_val.*phi);
        M_a_s_pmb = convectionUpwindTerm(-dlabda_w_face.*cp_face.*grad_p0, ut);
        % M_a_cp, RHS_a_cp = convectionTvdTerm(uw, cp_val, FL, ut)
        M_d_p_pmb = diffusionTerm(-labda_w_face.*cp_face);
        M_a_cs_pmb = convectionUpwindTerm(-k_face.*krw_face.*d_cp_mu_dcs.*grad_p0, ut);
        M_a_cp_pmb = convectionUpwindTerm(-k_face.*krw_face.*d_cp_mu_dcp.*grad_p0, ut);

        % create the PDE system M (p;s;c)=RHS
        % x_val = (p_val.value(:); sw_val.value(:); c_val.value(:))
        M = [M_bc_p+M_d_p_wmb,   M_t_s_wmb+M_a_s_wmb,  M_a_cs_wmb,  M_a_cp_wmb;
            M_d_p_omb,   M_bc_s+M_t_s_omb+M_a_s_omb,  zeros((Nx+2)*(Ny+2),(Nx+2)*(Ny+2)),  zeros((Nx+2)*(Ny+2),(Nx+2)*(Ny+2));
            M_d_p_smb,  M_a_s_smb+M_t_s_smb,  M_a_cs_smb+M_t_cs_smb+M_bc_cs,  M_a_cp_smb;
            M_d_p_pmb,  M_a_s_pmb+M_t_s_pmb,  M_a_cs_pmb,  M_t_cp_pmb+M_bc_cp+M_a_cp_pmb];

        RHS = [RHS_bc_p+RHS_t_s_wmb+M_a_s_wmb*sw_val.value(:)+M_a_cs_wmb*cs_val.value(:)+M_a_cp_wmb*cp_val.value(:); ...
               RHS_bc_s+RHS_t_s_omb+(M_a_s_omb)*sw_val.value(:); ...
               RHS_t_cs_smb+RHS_t_s_smb+RHS_bc_cs+M_a_s_smb*sw_val.value(:)+M_a_cs_smb*cs_val.value(:)+M_a_cp_smb*cp_val.value(:); ...
               RHS_t_cp_pmb+RHS_t_s_pmb+RHS_bc_cp+M_a_s_pmb*sw_val.value(:)+M_a_cs_pmb*cs_val.value(:)+M_a_cp_pmb*cp_val.value(:)];

        % x_sol = solveLinearPDE(m, M, RHS)
        
        try
            x_sol = M\RHS;
        catch
            loop_countor = max_int_loop;
            continue
        end

        p_new = reshape(x_sol(1:(Nx+2)*(Ny+2)), Nx+2, Ny+2);
        s_new = reshape(x_sol((1:(Nx+2)*(Ny+2))+(Nx+2)*(Ny+2)), Nx+2, Ny+2);
        cs_new = reshape(x_sol((1:(Nx+2)*(Ny+2))+2*(Nx+2)*(Ny+2)), Nx+2, Ny+2);
        cp_new = reshape(x_sol((1:(Nx+2)*(Ny+2))+3*(Nx+2)*(Ny+2)), Nx+2, Ny+2);

        error_s = sumabs(s_new(2:end-1, 2:end-1)-internalCells(sw_val));
        error_cp = sumabs(cp_new(2:end-1, 2:end-1)-internalCells(cp_val));

        %   error_c = sum(abs, internalCells(cp_val)-internalCells(cp_val))
        % println(error_s)
        % println(error_c)
        % w_p = 0.9
        % w_sw = 0.8
        % w_c = 0.8
        % p_val.value(:) = w_p*p_new(:)+(1-w_p)*p_val.value(:)
        % sw_val.value(:) = w_sw*s_new(:)+(1-w_sw)*sw_val.value(:)
        % c_val.value(:) = w_c*c_new(:)+(1-w_c)*c_val.value(:)
        p_val.value(:) = p_new(:);
        dsw_apple = 0.2;
        eps_apple = sqrt(eps); % 1e-5
        for i = 1:numel(s_new)
            if s_new(i)>=(1-sor)
                if sw_val.value(i)<(1-sor-eps_apple)
                    sw_val.value(i) = 1-sor-eps_apple;
                else
                    sw_val.value(i) = 1-sor;
                end
            elseif s_new(i)<=swc
                if sw_val.value(i)> swc+eps_apple
                    sw_val.value(i) = swc+eps_apple;
                else
                    sw_val.value(i) = swc;
                end
            elseif abs(s_new(i)-sw_val.value(i))>dsw_apple
                sw_val.value(i) = sw_val.value(i)+dsw_apple*sign(s_new(i)-sw_val.value(i));
            else
                sw_val.value(i) = s_new(i);
            end
        end
        sw_val = createCellVariable(m, internalCells(sw_val), BCs);
        cs_new(cs_new<0.0) = 0.0;
        cp_new(cp_new<0.0) = 0.0;
        
        cs_val.value(:) = cs_new(:);
        cp_val.value(:) = cp_new(:);

        cs_val = createCellVariable(m, cs_new, BCcs);
        cp_val = createCellVariable(m, cp_new, BCcp);

        % GR.plot(sw_val.value)
        % sw_val.value(:) = w_sw*s_new(:)+(1-w_sw)*sw_val.value(:)
    end
    if loop_countor<max_int_loop
        p_init = p_val;
        sw_init = sw_val;
        cp_init = cp_val;
        cs_init = cs_val;
        
        t = t+dt;
        dt = dt0;
        disp(t/t_final*100)
        visualizeCells(sw_init); drawnow;
    end
end
% JLD.save('results2D/all_data_heterogen_wf.jld', 'p', p_val, 'c', c_val,
% 'c_oil', c_oil, 'sw', sw_val, 'c_tracer', c_t_val)
% figure(1)
figure(1)
visualizeCells(sw_init)
% title('Sw')
figure(2)
visualizeCells(cs_init)
% title('c_DME')
figure(3)
visualizeCells(cp_init)
% title('c_tracer')
% savefig('profile.png')
% figure(2)
% plot(t_s, rec_fact)
% savefig('recovery.png')
% figure(3)
% plot(t_s, dp_hist)
% savefig('dp.png')
% PyPlot.pcolormesh(sw_val.value(2:end-1, 2:end-1)', cmap='YlGnBu', shading = 'gouraud', vmin=0, vmax=1)

  