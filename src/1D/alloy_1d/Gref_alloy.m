function [G_ref, A, B, T_diff] = Gref_alloy(a, z, numAtom, R, offset, pot_radial, ...
                                                        R0, pot_radial_ref)
% Solve 1d eigenvalue problem (-d^2/dx^2 + V(x)) u(x) = Eu(x) within MST 
% Generate the structural Green's function of reference, T_diff, and 
% two integral values A and B (see the manuscript for details).
% The integrals are indexed with the angular momentum and related to 
% the regular and irregular local radial solutions.
% 
% 14/Sep/2023

N = 1e2;
t_mat = zeros(2, numAtom);     % t-matrix
A = zeros(2, 1);
B = zeros(2, 1);

% the central atom is A by default
% atom A
[t_temp, A_temp, B_temp] = local_solver_1d(1, z, R, pot_radial, N);    
t_mat(1, 1) = t_temp;
A(1) = A_temp;
B(1) = B_temp;    
[t_temp, A_temp, B_temp] = local_solver_1d(-1, z, R, pot_radial, N);    
t_mat(2, 1) = t_temp;
A(2) = A_temp;
B(2) = B_temp;    

for j = 2 : numAtom
    t_mat(:, j) = t_mat(:, 1);
end

% atom B
tmat_B = zeros(2,1);
V = @(x) 2*pot_radial(x);
[t_temp, ~, ~] = local_solver_1d(1, z, R, V, N);    
tmat_B(1, 1) = t_temp; 
[t_temp, ~, ~] = local_solver_1d(-1, z, R, V, N);    
tmat_B(2, 1) = t_temp;

for j = 1 : (numAtom-1)/2
    idx = a(j);
    t_mat(:, idx) = tmat_B;
end
T = diag(reshape(t_mat, [ ], 1));

    
% The interpolation is used when z is on the real axis
if real(z) >0 && abs(imag(z))<=1e-2
    z1 = z + 1i*0.1;
    [G_ref1, T_ref1] = structural_green_1d(z1, numAtom, offset, pot_radial_ref, R0, N);    
    z2 = z - 1i*0.1;    
    [G_ref2, T_ref2] = structural_green_1d(z2, numAtom, offset, pot_radial_ref, R0, N);
    G_ref = (G_ref1+G_ref2)/2;
    T_ref = (T_ref1+T_ref2)/2;
else
    [G_ref, T_ref] = structural_green_1d(z, numAtom, offset, pot_radial_ref, R0, N);
end

% Check the decay rate of G_ref (all entries)
if (0)
x = zeros(numAtom, numAtom);
y = zeros(numAtom, numAtom);
for j = 1 : numAtom
    for n = 1 : numAtom
        x(j, n) = abs(j-n);
        y(j, n) = norm(G_ref(2*j-1:2*j, 2*n-1:2*n), 'fro');
    end
end
x = reshape(x, [ ], 1);
y = reshape(y, [ ], 1);
figure
semilogy(x, y, 'bo', 'markersize', 10)
end

T_diff = T - T_ref;
% M = T_diff^(-1) - G_ref;   % M-matrix used to be inverted


% Check the assumption \|G_ref T_diff\|<1
if (1)
Q = G_ref * T_diff;
Q_norm1 = norm(Q, 1);
fprintf('The 1-norm of contraction matrix is %f. \n', Q_norm1)

end

