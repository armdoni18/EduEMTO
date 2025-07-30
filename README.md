# 2D Magnetic Actuator Topology Optimization (Electromagnetics)

This MATLAB code performs 2D topology optimization of a magnetic actuator using a finite element method (FEM) with vector potential formulation. The objective is to maximize the magnetic force using Maxwell Stress Tensor (MST) under volume constraints and SIMP penalization.

---

## Features

- 2D FEM simulation using triangular mesh
- Gmsh-based geometry input (".msh")
- SIMP-based topology optimization with three-field projection
- Magnetic field solution with vector potential 
- Magnetic force computation via MST
- Adjoint sensitivity analysis
- Optimization with MMA (Method of Moving Asymptotes)
- Automatic plotting and saving of results

---

## Folder Structure

- Main_node_EM_TO.m                  % Main optimization script
- Magnetic_Actuator_Fine_Mesh.msh    % Gmsh mesh input (v4.1)
- Func2_pre_msh_2D.m                 % Mesh parser for Gmsh output
- Func3_Pre_FEM_init.m               % FEM system initialization
- Func4_Pre_Opt_init.m               % Optimization parameters and filter setup
- Func5_VecPot_A.m                   % FEM solver for vector potential A_z
- Func6_MagFlux_BxBy.m               % Compute magnetic flux density Bx, By, and |B|
- Func7_MST_Force.m                  % Compute magnetic force (Fy) using MST
- Func8_Sens.m                       % Sensitivity analysis (Adjoint + Heaviside projection)
- Func9_Opt_Plot.m                   % Plotting and result visualization


---

## How to Run

1. **Set the Mesh File in the Code:**

    In the **Main_node_EM_TO.m** file, locate and set the mesh filename in the following line:
    ```
    modelname = 'Magnetic_Actuator_xxxxxx';
    ```
    
    - Choose between the two mesh options:
        - `Magnetic_Actuator_Fine_Mesh.msh`           → full model
        - `Magnetic_Actuator_Fine_Mesh_Symmetry.msh`  → half model (symmetry)
    
    Example — to use the symmetry-based model:
    ```
    modelname = 'Magnetic_Actuator_Fine_Mesh_Symmetry';
    ```

2. **Adjust Volume Parameters (if needed):**

    - If using the full model, use:
    ```
    inputs.VT  = 97500;
    inputs.VND = 72150;
    inputs.VDD = 25350;
    ```
    
    - If using the half model, divide all three values by 2:
    ```
    inputs.VT  = 48750;
    inputs.VND = 36075;
    inputs.VDD = 12675;
    ```
  
3. **Update Boundary Conditions**
In the file **Func3_Pre_FEM_init.m**, the boundary conditions are set manually based on the domain size. You must adjust them based on whether you are using the full or symmetry mesh. If using the symmetry model, you must update this section to exclude the right boundary from Dirichlet conditions:

    - For the **full model** (Magnetic_Actuator_Fine_Mesh.msh), use:
    ```
    left_ind    = find(abs(fem.X(:,1)-(+0))<1e-6);
    right_ind   = find(abs(fem.X(:,1)-(+390))<1e-6);
    top_ind     = find(abs(fem.X(:,2)-(+0))<1e-6);
    bottom_ind  = find(abs(fem.X(:,2)-(+250))<1e-6);
     
    fem.bcdof   = unique([left_ind; right_ind; top_ind; bottom_ind]);
    fem.bcval   = zeros(length(fem.bcdof),1);
    ```

    - For the **symmetry model** (Magnetic_Actuator_Fine_Mesh_Symmetry.msh) use:
      - Exclude the right boundary from Dirichlet BCs:
    ```     
    left_ind    = find(abs(fem.X(:,1)-(+0))<1e-6);
    % right_ind   = find(abs(fem.X(:,1)-(+390))<1e-6); % Not used in symmetry
    top_ind     = find(abs(fem.X(:,2)-(+0))<1e-6);
    bottom_ind  = find(abs(fem.X(:,2)-(+250))<1e-6);
          
    fem.bcdof   = unique([left_ind; top_ind; bottom_ind]);
    fem.bcval   = zeros(length(fem.bcdof),1);
    ```
      - Treat the right boundary as a Neumann boundary with natural condition ∂𝐴𝑧/∂𝑛 = 𝑔 = 0

4. **Create Output Folder:**

    Users must create a folder named **Figures/** in the working directory.
    All iteration-based result plots will be saved there.
    
    **Important: The code will not run if the Figures/ folder is missing.**

5. **Run the Code in MATLAB:**

    Launch MATLAB, navigate to the project folder, and run:
    
    - Main_node_EM_TO.m
    
    The optimization will begin. Plots will be automatically saved in **Figures/** at each iteration.
