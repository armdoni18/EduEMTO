clear; clc; close all; tic;
modelname = 'Magnetic_Actuator_Fine_Mesh';
   
    %% Pre-Processing
    inputs.penal    = 3.0;          % SIMP penalization factor
    inputs.initdv   = -0.1;         % Initial design variable
    inputs.VT       = 97500;        % Total volume full model
    inputs.VND      = 72150;        % Non-design domain volume (Fixed) full model  
    inputs.VDD      = 25350;        % Design domain volume full model
    inputs.volfrac  = 0.25;         % Targeted volume fraction
    inputs.u0       = 4*pi*1e-7;    % Vacuum permeability
    inputs.mprop_D  = 1500;         % Reluctivity material in design domain (iron)
    inputs.mprop_A  = 1;            % Reluctivity of air
    inputs.mprop_C1 = 1;            % Reluctivity of coil 1 (+J)
    inputs.mprop_C2 = 1;            % Reluctivity of coil 2 (-J)
    inputs.mprop_ND = 1500;         % Reluctivity of C-core (iron)
    inputs.J_am2    = 1.5e2;        % Current density [A/m^2]
    inputs.conv     = 0.00008;      % Convergence criteria
    inputs.bt_init  = 0.1;          % Initial beta value for Heaviside projection
    inputs.bt_ic    = 1.5;          % Increment factor beta 
    inputs.bt_ns    = 4;            % Number of step for beta
    inputs.bt_fn    = 10;           % Final value for beta
    inputs.MMA_c    = 1000;         % MMA constant parameter
    inputs.rmin     = 2.0;            % Helmholtz filter radius              
    
    [mesh]   = Func2_pre_msh_2D(modelname);             % Load mesh information

    [fem]    = Func3_Pre_FEM_init(inputs,mesh);         % FEM set up

    [opt]    = Func4_Pre_Opt_init(inputs,fem);          % Optimization set up     

    fprintf('Pre-Processing Completed ✅\n');
    fprintf('Elapsed time for Pre-Processing: %s\n', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'))
       
    %% Main processing - Topology optimization
    while opt.bt < inputs.bt_fn

        % Three-field projection and element density calculation
        opt.fdv     = opt.chol_Kft'\(opt.chol_Kft\(opt.Tft*opt.nv)); 
        opt.nrho    = max(min(tanh(opt.bt*opt.fdv)/(2*tanh(opt.bt))+0.5,1),-1);
        opt.erho    = opt.Ten*opt.nrho;

        [fem]                                               = Func5_VecPot_A(fem, opt, inputs);                              
        [fem.Bx, fem.By, fem.B, fem.MagEn, fem.TotalMagEn]  = Func6_MagFlux_BxBy(fem);  
        [fem.Fy_total]                                      = Func7_MST_Force(fem);                                  
        [f, g, dfdx, dgdx, dfdrho_e]                        = Func8_Sens(fem, opt);                     

        % Store history
        opt.fhis(opt.iter) = f;                                                                      
        opt.ghis(opt.iter) = g;                                                                      

        % Convergence check
        if (opt.iter > 1)
            opt.deltaf = abs((opt.fhis(opt.iter) - opt.fhis(opt.iter - 1)) / opt.fhis(opt.iter - 1));
        end    

        % Display & plot
        disp(fprintf('Iter:%d, f:%.4f, Volume:%.4f, deltaf:%.5f, beta:%.2f', ...
             opt.iter, f, g + opt.volfrac, opt.deltaf, opt.bt));
        Func9_Opt_Plot(fem,opt); 

        % MMA update
        [dvnew,~,~,~,~,~,~,~,~,opt.MMA.low,opt.MMA.upp] = mmasub(...
            1,length(opt.dv),opt.iter,opt.dv,opt.dvmin,opt.dvmax,opt.dvold,opt.dvolder,...
            f,dfdx,g,dgdx,opt.MMA.low,opt.MMA.upp,opt.MMA.a0,opt.MMA.a,opt.MMA.c,opt.MMA.d);        

        opt.iter             = opt.iter + 1;       
        opt.dvolder          = opt.dvold; 
        opt.dvold            = opt.dv; 
        opt.dv               = dvnew;      
        opt.nv(opt.dof_dd,1) = opt.dv;
    
        % Continuation
        if (exist('cont_sw','var') == 0) && (opt.deltaf < inputs.conv)
            cont_sw = 1; cont_iter = 0;
            disp('Continuation for Heaviside Projection')
        elseif exist('cont_sw','var')
            cont_iter = cont_iter + 1;
            if mod(cont_iter, inputs.bt_ns) == 1
                opt.bt = opt.bt * inputs.bt_ic;
            end
        end    
    end

    fprintf('Main-Processing Completed ✅\n');
    fprintf('Elapsed time for Main Processing:%s\n',(datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')))
 %%   
    saveas(figure(1), 'Figures/Obj_Final.png');                                                 % Save the final objective function plot
    saveas(figure(2), 'Figures/Vol_Final.png');                                                 % Save the final volume constraint plot
% %% Post-processing
%     Func10_PostProcessold(fem,opt,inputs); 
%     fprintf('Elapsed time for Post-Processing:%s\n',(datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')))
%     fprintf('Post-Processing Completed ✅\n');