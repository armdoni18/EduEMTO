function [mesh] = Func2_pre_msh_2D(modelname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective
%   Convert "filename.msh(v4.1)" to "Inner FEM data(IX,X,bc,lc)"
%
% Requirement
%   (GMSH) Physical tag is required to identify domain(Design, Non-design,Coil, and Air)
%
%  ***Domain identification***
%
%  If there are no non-design domain, any physical tag is not required to 
%  identify design domain. Because default domain setting is design domain.
%  However if there are a non-design domian, physical tag is required. A
%  physical tag must contain specific keyword for target domain.
%
%  'NonDesign' is keywords for check non-design domain 
%                                  
%  ***************************
%  
% IX : element conectivity + Domain index information
%  Type : index matrix
%    IX(:,1) : 1st node index
%    IX(:,2) : 2nd node index
%    IX(:,3) : 3rd node index
%    IX(:,4) : Domain index      % 1 : Air // 2 : Design
%                                % 3 : Coil1  // 4 : Coil2   // 5 : NonDesign
%  
% X : Node coordinates information
%  Type : coordinates matrix
%    X(:,1) : X-axis coordinated
%    X(:,2) : Y-axis coordinated
%
% lc : Loading condition {Cell}(element conectivity matrix)
%   lc.c{ind} : Current loading
%   lc.b{ind} : Body loading
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fid = fopen([modelname,'.msh']);
flag_category = 0; % Flag for ctaegroy
flag_state = 0;    % Flag for state

while 1
    % File reading part
    tline = fgetl(fid);
    if ~ischar(tline), break, end % if end the file, break the while loop
    
    % line 
%     [data, tf] = str2num(tline);
    [data,~,errmsg,nextindex] = sscanf(tline,'%f');
%     disp(tline)
%     disp(data)
    
    if (~isempty(errmsg)) && (nextindex == 1)
        if strcmp(tline,'$PhysicalNames')
          flag_category = 1;
        elseif strcmp(tline,'$Entities')
          flag_category = 2;
        elseif strcmp(tline,'$Nodes')
          flag_category = 4;
        elseif strcmp(tline,'$Elements')
          flag_category = 5;
        end
        flag_state = 0;
        continue
    end
    
    if flag_category == 1 % $PhysicalNames
        if flag_state == 0
            if length(data) ~= 1
                error
            end
            num_phys_entities = data;
            flag_state = 1;
%%%%%% Check %%%%%%
            phys_tag_str.Air = [];
            phys_tag_str.Design = [];
            phys_tag_str.Coil1 = [];
            phys_tag_str.Coil2 = [];
            phys_tag_str.NonDesign = [];
%%%%%%%%%%%%%%%%%%%
            
            ind_phys_entities = 0;
        
                                                                           
        elseif flag_state == 1  
%%%%%% Check %%%%%%
            if (length(strfind(tline(nextindex:end),'"Air"')) == 1)
                Air_list = [];
                ind_phys_entities = ind_phys_entities + 1;
                phys_tag_str.Air = [phys_tag_str.Air, (data(2))];
            elseif (length(strfind(tline(nextindex:end),'"Design"')) == 1)
                Design_list = [];
                ind_phys_entities = ind_phys_entities + 1;
                phys_tag_str.Design = [phys_tag_str.Design, (data(2))];
            elseif (length(strfind(tline(nextindex:end),'"Coil1"')) == 1)
                Coil1_list = [];
                ind_phys_entities = ind_phys_entities + 1;
                phys_tag_str.Coil1 = [phys_tag_str.Coil1, (data(2))];
            elseif (length(strfind(tline(nextindex:end),'"Coil2"')) == 1)
                Coil2_list = [];
                ind_phys_entities = ind_phys_entities + 1;
                phys_tag_str.Coil2 = [phys_tag_str.Coil2, (data(2))];
            elseif (length(strfind(tline(nextindex:end),'"NonDesign"')) == 1)
                NonDesign_list = [];
                ind_phys_entities = ind_phys_entities + 1;
                phys_tag_str.NonDesign = [phys_tag_str.NonDesign, (data(2))];
            end
%%%%%%%%%%%%%%%%%%%
        end
        
        if ind_phys_entities >= num_phys_entities
            flag_state = 0;
            ind_phys_entities = 0;
        end
                    
    elseif flag_category == 2 % $Entities
        if flag_state == 0
            num_point = data(1);
            num_line = data(2);
            num_surf = data(3);
            num_vol = data(4);
            if exist('phys_tag_str','var')
                if ~(isempty(phys_tag_str))
                    flag_state = 1;
                    flag_entity = 0;
                    ind_entity = 0;
                end
            end
            
        elseif flag_state == 1
            if flag_entity == 0 % Point
                ind_entity = ind_entity + 1;
                if ind_entity >= num_point
                    flag_entity = 1;
                    ind_entity = 0;
                end
            end
            
            if flag_entity == 1 % Curve
                ind_entity = ind_entity + 1;
                                
                if ind_entity >= num_line
                    flag_entity = 2;
                    ind_entity = 0;
                end
            end
            
            if flag_entity == 2 % Surface
                
                ind_entity = ind_entity + 1;
%%%%%% Check %%%%%%
                for i = 1:data(8)
                    if (sum(data(8+i) == phys_tag_str.Air) >= 1)                % Air
                        Air_list = [Air_list; data(1)];
                    elseif (sum(data(8+i) == phys_tag_str.Design) >= 1)    	    % Design
                        Design_list = [Design_list; data(1)];
                    elseif (sum(data(8+i) == phys_tag_str.Coil1) >= 1)          % Coil1
                        Coil1_list = [Coil1_list; data(1)];
                    elseif (sum(data(8+i) == phys_tag_str.Coil2) >= 1)          % Coil2
                        Coil2_list = [Coil2_list; data(1)];
                    elseif (sum(data(8+i) == phys_tag_str.NonDesign) >= 1)      % Nondesign
                        NonDesign_list = [NonDesign_list; data(1)];
                    end
                end                          
%%%%%%%%%%%%%%%%%%%          
                
                if ind_entity >= num_surf
                    flag_entity = 3;
                    ind_entity = 0;
                end
            end
            
            if flag_entity == 3 % Volume
                
                ind_entity = ind_entity + 1;
                
                if ind_entity >= num_vol
                    flag_entity = 4;
                    ind_entity = 0;
                end
            end
        end
        
    elseif flag_category == 4           % $Nodes
        if flag_state == 0              % Initial state
            num_node = data(2);         % Total number of node
            ind_node = 0;               % Node index initializing
            p = zeros(num_node,2);      % Node coordinates data
            flag_state = 1; %
        
        elseif flag_state == 1          % Info data parsing
            num_loop = data(4);         % Number of node on this block
            if num_loop == 0
                continue
            else
                flag_state = 2;
                ind_inner_loop = 0;
            end
            ind_list = zeros(num_loop,1);
            
        elseif flag_state == 2                      % ind data parsing
            ind_inner_loop = ind_inner_loop + 1;
            
            ind_list(ind_inner_loop) = data;
            
            if ind_inner_loop >= num_loop           % If ind data parsing is ended
                flag_state = 3;
                ind_inner_loop = 0;
            end
            
        elseif flag_state == 3                      % point data parsing
%            ind_node = ind_node + 1;
            ind_inner_loop = ind_inner_loop + 1;
            try
                p(ind_list(ind_inner_loop),:) = data(1:2);
            catch
                disp(ind_node)
            end
            
            if ind_inner_loop >= num_loop
                flag_state = 1;
                ind_inner_loop = 0;
            end
        elseif flag_state == 99                 % error exception
            ind_inner_loop = ind_inner_loop + 1;
            
            if ind_inner_loop >= num_loop*2     % If ind data parsing is ended
                flag_state = 1;
                ind_inner_loop = 0;
            end
        end
    elseif flag_category == 5           % $Elements
        if flag_state == 0              % Initial state
            if length(data) ~= 4
                error
            end
            ind_elem = 0; 
            flag_state = 1;     %
            t = [];             % Element conectivity information
            
        elseif flag_state == 1  % Info data parsing
            num_loop = data(4);    % Number of node on this block
            flag_state = 99;
            ind_inner_loop = 0;
            
            flag_condition = 0;
            
            if data(1) == 1      % 1D element : Line
                continue
            elseif data(1) == 2   % 2D element
                if data(3) == 2  % 2D CST
                    flag_state = 4;
                    num_elem = data(4);
                    t_temp = zeros(num_elem,4);
%%%%%% Check %%%%%%
                    if exist('Air_list','var')                       
                        if sum(Air_list==data(2))
                            t_temp(:,4) = ones(num_elem,1);
                        end
                    end
                    
                    if exist('Design_list','var')                       
                        if sum(Design_list==data(2))
                            t_temp(:,4) = ones(num_elem,1)*2;
                        end
                    end
                    
                    if exist('Coil1_list','var')                       
                        if sum(Coil1_list==data(2))
                            t_temp(:,4) = ones(num_elem,1)*3;
                        end
                    end

                    if exist('Coil2_list','var')                       
                        if sum(Coil2_list==data(2))
                            t_temp(:,4) = ones(num_elem,1)*4;
                        end
                    end
                    
                    if exist('NonDesign_list','var')                       
                        if sum(NonDesign_list==data(2))
                            t_temp(:,4) = ones(num_elem,1)*5;
                        end
                    end
%%%%%%%%%%%%%%%%%%%

                elseif data(3) == 3 % 2D quadrangle element
                    continue
                end
            end
            
        elseif flag_state == 2  % 1D boundary element conectivity data parsing
            ind_inner_loop = ind_inner_loop + 1;
            
                       
            if ind_inner_loop >= num_loop
                flag_state = 1;
                ind_inner_loop = 0;
            end
            
        elseif flag_state == 3  % 1D boundary element conectivity data parsing
            ind_inner_loop = ind_inner_loop + 1;
            
                       
            if ind_inner_loop >= num_loop
                flag_state = 1;
                ind_inner_loop = 0;
            end
            
        elseif flag_state == 4  % 2D CST element conectivity data parsing
            ind_inner_loop = ind_inner_loop + 1;
            t_temp(ind_inner_loop,1:3) = data(2:4);
            
            if ind_inner_loop >= num_loop
                flag_state = 1;
                ind_inner_loop = 0;
                t = [t; t_temp];
            end
            
        elseif flag_state == 99 % error exception
            ind_inner_loop = ind_inner_loop + 1;
            
            if ind_inner_loop >= num_loop % If ind data parsing is ended
                flag_state = 1;
                ind_inner_loop = 0;
            end
        end
    end
end
fclose(fid);
%% Output setting
mesh.IX = t;
mesh.X = p;

fprintf('Mesh Generation Done. ✅\n');
end