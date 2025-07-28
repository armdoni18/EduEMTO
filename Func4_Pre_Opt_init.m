function [opt,MMA] = Func4_Pre_Opt_init(inputs,fem)
% Separation of non-design domains
    values = [1, 3, 4, 5];  
    nd_ele = find(ismember(fem.IX(:,4), values));
    opt.dof_nd      = unique([fem.IX(nd_ele,1);fem.IX(nd_ele,2);fem.IX(nd_ele,3)]); % DOF in non-design domain
    opt.dof_dd      = setdiff((1:fem.nn)',opt.dof_nd);                              % DOF in desing domain 
    opt.nn_dd       = size(opt.dof_dd,1);                                           % Number of nodes in design domain
% Separation of design domains
    dd_ele          = 2;                                                            % Design domain element
    opt.erho_dd     = find(ismember(fem.IX(:,4), dd_ele));                          % List of element in design domain
    opt.ne_dd       = size(opt.erho_dd,1);                                          % Number of element in design domain   
% Setting problem parameters
    opt.VT = inputs.VT;                                                             % Total volume (include all components)
    opt.VND = inputs.VND;                                                           % Non-design domain volume
    opt.VDD = inputs.VDD;                                                           % Design domain volume design
    opt.volfrac = inputs.volfrac;                                                   % Volume fraction setting
    opt.penal = inputs.penal;                                                       % Penalization parameter   
% Setting of optimization variables    
    opt.dv = ones(length(opt.dof_dd),1) * inputs.initdv;                            % Design domain node density (influced by inital value)
    opt.nv(opt.dof_nd,1) = 1;                                                       % Non-deisgn domain node density value (fixed)
    opt.nv(opt.dof_dd,1) = opt.dv;
    opt.dvold=opt.dv;           opt.dvolder=opt.dv;            
    opt.dvmin=opt.dv*0-1;        opt.dvmax=opt.dv*0+1;
    opt.iter=1;                   
    opt.deltaf = 1.0;
    opt.deltaf2 = 1.0;
    opt.deltaf3 = 1.0;  
    opt.bt = inputs.bt_init ;               
% Setting of MMA variables
    MMA.a0 = 1;                 MMA.a = 0;   
    MMA.c = [inputs.MMA_c;];    MMA.d = 1;
    MMA.low = opt.dvmin;        MMA.upp = opt.dvmax; 
    opt.MMA = MMA;
%% Build filter element stiffness (Kft) and Transformation matrix (Tft) 
    nx = [fem.X(fem.IX(:,1),1) fem.X(fem.IX(:,2),1) fem.X(fem.IX(:,3),1)];          % Element node x location     
    ny = [fem.X(fem.IX(:,1),2) fem.X(fem.IX(:,2),2) fem.X(fem.IX(:,3),2)];          % Element node y location     
    Kd = (inputs.rmin/2/sqrt(2))^2 * [1 0; 0 1];                              
    NN = [2 1 1; 1 2 1; 1 1 2]/6;                                               
    isf = reshape((kron(fem.IX(:,1:3),ones(1,3)))',9*fem.ne,1)';
    jsf = reshape((kron(fem.IX(:,1:3),ones(3,1)))',9*fem.ne,1)';
    for e = 1:fem.ne
        px = [nx(e,1) nx(e,2) nx(e,3)]';        
        py = [ny(e,1) ny(e,2) ny(e,3)]';        
        Ve = fem.Ve(e);
        beta_i  = py(2) - py(3);
        beta_j  = py(3) - py(1);
        beta_m  = py(1) - py(2);
        gamma_i = px(3) - px(2);
        gamma_j = px(1) - px(3);
        gamma_m = px(2) - px(1);
        B = 1/(2*Ve) * [beta_i, beta_j, beta_m;...
                       gamma_i, gamma_j, gamma_m];
        Kft((e-1)*9+1:(e)*9) = reshape((B'*Kd*B + NN)*Ve,9,1);                  
        Tft((e-1)*9+1:(e)*9) = reshape(NN * Ve, 9, 1);                          
    end
opt.chol_Kft = chol(sparse(isf, jsf, Kft), 'lower');
opt.Tft = sparse(isf, jsf, Tft);
%% Build matrix for transformation from nodal to element density (Ten)
    opt.Ten=sparse(fem.ne,fem.nn);
    for e=1:fem.ne
       opt.Ten(e,fem.IX(e,1:3))=1/3;
    end
fprintf('Optimization Initialization Done. ✅\n');
end