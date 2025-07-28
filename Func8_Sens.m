function [f, g, dfdx, dgdx, dfdrho_e] = Func8_Sens(fem,opt)
%% Compute sensitivity of magnetic force Fy w.r.t. rho
f = fem.Fy_total;
%% === Step 1: Compute dF/dA ===
dfdA = zeros(fem.ndof,1);
IntegPathEdges  = size(fem.cleaned_air_loop_around_plunger, 1);
plunger_nodes   = unique(fem.IX(fem.IX(:,4)==2, 1:3));
plunger_coords  = fem.X(plunger_nodes, :);
plunger_center  = mean(plunger_coords, 1);

for k = 1:IntegPathEdges
    n1 = fem.cleaned_air_loop_around_plunger(k,1);
    n2 = fem.cleaned_air_loop_around_plunger(k,2);
    p1 = fem.X(n1, :);
    p2 = fem.X(n2, :);
    seg = p2 - p1;
    ds  = norm(seg);
    normal = [seg(2), -seg(1)];
    normal = normal / norm(normal);

    mid = (p1+p2)/2;

    if dot(normal, plunger_center - mid) > 0
        normal = -normal;
    end

    % Shift midpoint slightly inward (toward plunger)
    eps_shift = 1e-3;
    shifted_mid = mid + eps_shift * normal;

    % Find element whose centroid is closest to the shifted point
    min_dist = inf;     closest_e = 1;
    for e = 1:fem.ne
        nodes_e = fem.IX(e, 1:3);
        coords = fem.X(nodes_e, :);
        centroid = mean(coords, 1);
        dist = norm(shifted_mid - centroid);
        if dist < min_dist
            min_dist = dist;
            closest_e = e;
        end
    end

    nodes_e = fem.IX(closest_e,1:3);
    coords = fem.X(nodes_e,:);
    Ve = fem.Ve(closest_e);

    xi=coords(:,1); yi=coords(:,2);
    bi=[yi(2)-yi(3); yi(3)-yi(1); yi(1)-yi(2)];
    ci=[xi(3)-xi(2); xi(1)-xi(3); xi(2)-xi(1)];

    Bx = fem.Bx(closest_e);
    By = fem.By(closest_e);
    mu = fem.IX(closest_e,5);

    dBx_dA = (1/(2*Ve)) * ci;
    dBy_dA = (-1/(2*Ve)) * bi;

    dTxy_dA = mu*(By*dBx_dA + Bx*dBy_dA);
    dTyy_dA = mu*0.5*(2*By*dBy_dA - 2*Bx*dBx_dA);

    dFy_edge = (dTxy_dA*normal(1) + dTyy_dA*normal(2)) * ds;

    dfdA(nodes_e) = dfdA(nodes_e) + dFy_edge;
end

%% === Step 2: Solve adjoint problem ===
lambda = fem.S \ dfdA;

%% === Step 3: Compute dF/drho per element using adjoint ===
dfdrho_e = zeros(fem.ne,1);
penal = opt.penal;
V0 = fem.V0;
Vmin = fem.Vmin;

for e=1:fem.ne
    rhoe = opt.erho(e);
    penal_term = penal*(V0-Vmin)*rhoe^(penal-1);
    Se = reshape(fem.S_S((e-1)*9+1:e*9),3,3);
    Ae_vec = fem.A(fem.IX(e,1:3));
    lambda_vec = lambda(fem.IX(e,1:3));
    dfdrho_e(e) = lambda_vec'*(penal_term/V0)*Se*Ae_vec;
end
    
%% === Step 4: Project to nodal design variables ===
DerhoDnrho = opt.Ten';      %%
DnrhoDfdv  = (sech(opt.bt*opt.fdv)).^2 * opt.bt/(2*tanh(opt.bt));  

dfdx_all = DerhoDnrho * dfdrho_e;

dfdx = dfdx_all .* DnrhoDfdv;
dfdx = dfdx(opt.dof_dd);

%% Compute sensitivity of Volume constraint
%% === Step 1: Constraint value
g           = ((fem.Ve*opt.erho)-opt.VND)/(opt.VT-opt.VND)-opt.volfrac;

%% === Step 2: Compute dV/derho (element-based)
dVdrho_e    = fem.Ve'/opt.VT;                          

%% === Step 3: Chain rule - d(rho_e)/d(rho_nodal)
DerhoDnrho  = opt.Ten';                              

%% === Step 4: Node-wise sensitivity (before projection)
dgdx_all_nodes = DerhoDnrho * dVdrho_e;

%% === Step 5: Chain with Heaviside filter sensitivity
DnrhoDfdv   = (sech(opt.bt*opt.fdv)).^2*opt.bt/(2*tanh(opt.bt));
dgdx        = dgdx_all_nodes .* DnrhoDfdv;
dgdx        = dgdx(opt.dof_dd);                                             % only designable dofs
dgdx        = dgdx';

fprintf('Sensitivity Computation Done.✅\n');
end
