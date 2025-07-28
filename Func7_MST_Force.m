function [Fy_total] = Func7_MST_Force(fem)
%% Maxwell Stress Tensor Force Computation on design domain boundary
    Fy_total = 0;
    IntegPathEdges  = size(fem.cleaned_air_loop_around_plunger, 1);
    % Compute centroid of plunger region
    plunger_nodes   = unique(fem.IX(fem.IX(:,4)==2, 1:3));
    plunger_coords  = fem.X(plunger_nodes, :);
    plunger_center  = mean(plunger_coords, 1);                              % [x, y]
    % calculate the edge length
    for k = 1:IntegPathEdges
        n1 = fem.cleaned_air_loop_around_plunger(k,1);
        n2 = fem.cleaned_air_loop_around_plunger(k,2);
        p1 = fem.X(n1, :);
        p2 = fem.X(n2, :);
        seg = p2 - p1;
        ds  = norm(seg);
  
        % Compute the outward normal vector
        normal = [seg(2), -seg(1)];                                         % perpendicular to seg
        normal = normal / norm(normal);                                     % normalize
        mid = (p1 + p2) / 2;                                                % Midpoint of the edge  
        to_plunger = plunger_center - mid;                                  % Vector from midpoint to plunger center

        if dot(normal, to_plunger) > 0                                      % Check direction: if normal points toward plunger, flip it
            normal = -normal;
        end

        eps_shift = 1e-3;                                                   % Shift midpoint slightly inward (toward plunger)
        shifted_mid = mid + eps_shift * normal;

        % Find element whose centroid is closest to the shifted point
        min_dist = inf;        closest_e = 1;
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

        % Magnetic field components and permeability from closest element
        Bx = fem.Bx(closest_e);
        By = fem.By(closest_e);
        mu = fem.IX(closest_e, 5);
  
        % Maxwell stress tensor
        T = (mu) * [0.5*(Bx^2 - By^2), Bx*By;
                    Bx*By, 0.5*(By^2 - Bx^2)];

        fem.T_cleaned_air_edges(:,:,k) = T;
    
        dF = (T * normal') * ds;                                            % Force contribution from this edge
   
        Fy_total = Fy_total + dF(2);                                        % Accumulate total force

    end
%%
    fem.Fy_total    = Fy_total;                                             % Objective function = vertical force
    
        fprintf('Magnetic force by MST computation Done.✅\n');
end