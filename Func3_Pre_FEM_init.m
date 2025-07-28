function [fem]    = Func3_Pre_FEM_init(inputs,mesh)
%% Setting up the FEM variables
    fem.IX      = mesh.IX;                                                      % Node indices and domain index 
    fem.X       = mesh.X;                                                       % Node coordinates
    fem.nn      = size(fem.X,1);                                                % Number of nodes amount
    fem.ndof    = fem.nn;                                                       % Number of degrees of freedom; A in Z direction only
    fem.ne      = size(fem.IX,1);                                               % Number of elements amount
    fem.edof    = [fem.IX(:,1) fem.IX(:,2) fem.IX(:,3)];                        % DOF information of element    
    fem.is      = reshape((kron(fem.edof,ones(1,3)))',9*fem.ne,1)';             % Row index vector of element stiffnes matrix
    fem.js      = reshape((kron(fem.edof,ones(3,1)))',9*fem.ne,1)';             % Coloumn index vector of element stiffness matrix
    fem.V0      = 1 / (inputs.u0 * inputs.mprop_D(1));                          % Solid (iron) reluctivity 
    fem.Vmin    = 1 / (inputs.u0 * inputs.mprop_A(1));                          % Void (air) reluctivity
    fem.D       = eye(2);                                                       % Identity matrix for element stiffness matrix formulation
%% Building the solid reluctivity matrix (S_S) and volume vector () 
    fem.nx=[fem.X(fem.IX(:,1),1) fem.X(fem.IX(:,2),1) fem.X(fem.IX(:,3),1)];    % x coordinate node each element
    fem.ny=[fem.X(fem.IX(:,1),2) fem.X(fem.IX(:,2),2) fem.X(fem.IX(:,3),2)];    % y coordinate node each element  
    for e=1:fem.ne
        px=[fem.nx(e,1) fem.nx(e,2) fem.nx(e,3)]';                              % Triangle x coodinate
        py=[fem.ny(e,1) fem.ny(e,2) fem.ny(e,3)]';                              % Triangle y coordinate
        A = abs(1/2*det([ones(3,1), px, py]));                                  % Triangular area computation
        % Shape function gradient    
        beta_i  = py(2)-py(3);
        beta_j  = py(3)-py(1);
        beta_m  = py(1)-py(2);
        gamma_i = px(3)-px(2);
        gamma_j = px(1)-px(3);
        gamma_m = px(2)-px(1);
        B = 1/(2*A)* [ beta_i, beta_j,beta_m;...                                % Matrix gradient  
                       gamma_i,gamma_j,gamma_m];
        fem.S_S((e-1)*9+1:(e)*9) = reshape(B'*fem.D*B*A,9,1);                   % Element stifness matrix 
        fem.Ve(e)= A;                                                           % storing the triangle element area
    end
%% Setting of FEM boundary condition index
    left_ind    = find(abs(fem.X(:,1)-(+0))<1e-6);                              % left boundary setting 
    right_ind   = find(abs(fem.X(:,1)-(+390))<1e-6);                            % right boundary setting
    top_ind     = find(abs(fem.X(:,2)-(+0))<1e-6);                              % top boundary setting         
    bottom_ind  = find(abs(fem.X(:,2)-(+250))<1e-6);                            % buttom boundary setting
    fem.bcdof   = unique([left_ind; right_ind; top_ind; bottom_ind]);           % BC node          
    fem.bcval   = zeros(length(fem.bcdof),1);                                   % BC node value = 0
%% Setting up the current in coil-1 domain
    fem.ncoil1 = [];
    for e = 1:size(mesh.IX, 1)
        if mesh.IX(e,4) == 3                                                    % check if domain is coil-1 
            coil1nodes = mesh.IX(e,1:3);                                        % get the triangular nodes for coil-1 domain
            fem.ncoil1 = [fem.ncoil1; coil1nodes];                              % store the triangular nodes for coil-1 domain
            mesh.ncoil1 = fem.ncoil1;  
        end
    end
% Setting of value current in coil-1 domain 
    fem.c_A_all = zeros(size(mesh.ncoil1, 1), 1);  
    fem.Tdof = [];      
    fem.Tval = [];  
    for j=1:size(mesh.ncoil1,1)
        v1 = fem.X(mesh.ncoil1(j,1), 1:2);                                      % get the triangle vertices for coil-1 domain
        v2 = fem.X(mesh.ncoil1(j,2), 1:2);
        v3 = fem.X(mesh.ncoil1(j,3), 1:2); 
        c_A = abs(det([v1 - v3; v2 - v3])) / 2;                                 % compute the area of the triangle in coil-1 domain
        fem.c_A_all(j) = c_A;                                                   % storing every triangle area in coil-1 domain
        % Assign current density to each vertex in coil-1
        for k=1:3
            fem.Tdof = [fem.Tdof; mesh.ncoil1(j,k)];
            fem.Tval = [fem.Tval; -c_A/3*inputs.J_am2];
        end
    end
%% Setting up the current in coil-2 domain
    fem.ncoil2 = [];  
    for e = 1:size(mesh.IX, 1)
        if mesh.IX(e,4) == 4                                                    % checked if domain is coil-2
            coil2nodes = mesh.IX(e,1:3);                                        % get the triangular nodes for coil-2 domain
            fem.ncoil2 = [fem.ncoil2; coil2nodes];                              % store the triangular nodes for coil-2 domain
            mesh.ncoil2 = fem.ncoil2;   
        end
    end
% Setting of value current in coil-2 domain 
    for j=1:size(mesh.ncoil2,1)
        v1 = fem.X(mesh.ncoil2(j,1), 1:2);                                      % get the triangle vertices for coil-2
        v2 = fem.X(mesh.ncoil2(j,2), 1:2);
        v3 = fem.X(mesh.ncoil2(j,3), 1:2); 
        c_A = abs(det([v1 - v3; v2 - v3])) / 2;                                 % compute the area of the triangle in coil-2 domain
        fem.c_A_all(j) = c_A;                                                   % storing every triangle area in coil-2 domain
        % Assign current density to each vertex in coil-2
        for k=1:3
            fem.Tdof = [fem.Tdof; mesh.ncoil2(j,k)];
            fem.Tval = [fem.Tval; c_A/3*inputs.J_am2];
        end
    end
%% Extract boundary segments for integration path for calculating the magnetic force by MST
    % === STEP 1: Get plunger boundary edges ===
    target_mat  = 2;                                                            % targeted material (design domain)  
    plunger_elements    = fem.IX(fem.IX(:,4) == target_mat, 1:3);               % desing domain element is detemined by number 2 (plunger)
    plunger_edges = [plunger_elements(:,[1 2]);
             plunger_elements(:,[2 3]);
             plunger_elements(:,[3 1])];
    plunger_edges = sort(plunger_edges, 2);                                     % for consistent ordering
    [unique_edges, ~, ic] = unique(plunger_edges, 'rows');                      % count occurrences
    counts = accumarray(ic, 1);
    fem.plunger_boundary_edges = unique_edges(counts == 1,:);                   % keep the outer plunger edges only
    % === STEP 2: Find all AIR triangles ===
    air_triangles = fem.IX(fem.IX(:,4) == 1, 1:3);                              % air triangles element
    % === STEP 3: Check if each air triangle touches the plunger boundary ===
    connected_air = [];
    for i = 1:size(air_triangles,1)
        tri_nodes = air_triangles(i,:);
        for j = 1:size(fem.plunger_boundary_edges,1)
            edge_nodes = fem.plunger_boundary_edges(j,:);
            if any(ismember(edge_nodes, tri_nodes))                             % if any node shared
                connected_air = [connected_air; tri_nodes];
                break;
            end
        end
    end
    % === STEP 4: Build edges from connected air triangles ===
    air_edges = [connected_air(:,[1 2]);
                 connected_air(:,[2 3]);
                 connected_air(:,[3 1])];
    air_edges = sort(air_edges, 2);                                             % consistent ordering
    [unique_aedges, ~, ic_air] = unique(air_edges, 'rows');
    aedge_counts = accumarray(ic_air, 1);
    outer_air_loop = unique_aedges(aedge_counts == 1, :);                       % boundary of air shell
    % === OUTPUT ===
    fem.air_loop_around_plunger = outer_air_loop;
    % Ensure both edge lists are sorted (should already be, but just in case)
    plunger_edges = sort(fem.plunger_boundary_edges, 2);
    air_loop_edges = sort(fem.air_loop_around_plunger, 2);
    % Use setdiff to remove overlapping edges
    fem.cleaned_air_loop_around_plunger = setdiff(air_loop_edges, plunger_edges, 'rows');
    fprintf('FEM Initialization Done. ✅\n');
end