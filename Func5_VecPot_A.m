function [fem] = Func5_VecPot_A(fem, opt, inputs)

% === 1. Initialize global system ===
S = sparse(fem.ndof, fem.ndof);
T = zeros(fem.ndof, 1);

% === 2. Initialize element-wise reluctivity vector ===
nu_e_all = zeros(fem.ne,1);

for e = 1:fem.ne
    domain_idx = fem.IX(e,4);

    if domain_idx == 2                                                      % Design domain
        rhoe = opt.erho(e);
        nu_e = fem.Vmin + (fem.V0 - fem.Vmin) * rhoe^opt.penal;
  
    elseif domain_idx == 5                                                  % Non-design domain
        nu_e = 1 / (inputs.u0 * inputs.mprop_ND(1));

    elseif domain_idx == 3                                                  % Coil-1
        nu_e = 1 / (inputs.u0 * inputs.mprop_C1(1));

    elseif domain_idx == 4                                                  % Coil-2
        nu_e = 1 / (inputs.u0 * inputs.mprop_C2(1));

    else                                                                    % Air
        nu_e = 1 / (inputs.u0 * inputs.mprop_A(1));
    end
    nu_e_all(e) = nu_e;
    fem.IX(e,5) = nu_e;

end

%% === 3. Build element-wise stiffness vector (Velist) ===
Velist = kron(nu_e_all, ones(9,1))' .* fem.S_S;

% === 4. Assemble global stiffness matrix ===
S = sparse(fem.is, fem.js, Velist);
S = (S + S') / 2;                                                           % Ensure symmetry for the matrix


% === 5. Apply boundary conditions ===
S(fem.bcdof, :) = 0;
S(:, fem.bcdof) = 0;
for i = 1:length(fem.bcdof)
    S(fem.bcdof(i), fem.bcdof(i)) = 1;
end
T(fem.bcdof) = fem.bcval;

% === 6. Apply current source ===
for i = 1:length(fem.Tdof)
    T(fem.Tdof(i)) = T(fem.Tdof(i)) + fem.Tval(i);
end

% === 7. Solve for vector potential A ===
A = S \ T;

fem.A = A;
fem.S = S;

fprintf('FEM Solved and Vector potential computation Done. ✅\n');
end
