function [G_ref, A, B, T_diff] = Gref_basis_chain(z, lmax, numAtom, R, offset, ...
    pot_radial, pot_radial_ref)
% Consider (-\nabla + V(r)) u(r) = Eu(r) within MST 
% Generate the structural Green's function of reference,
% two integral values A and B, and T_diff (see the manuscript for details), 
% which are indexed with the angular momentum.
% Two integral are related to the regular and irregular local radial solutions 
% for a muffin-tin system with identical potentials
%
% Note that the approximation of cell by muffin-tin ball is used in the
% computation of the integrals.
% 
% 09/Aug/2023

N = 1e2;
t_mat = zeros((lmax+1)^2, numAtom);     % t-matrix
A = zeros((lmax+1)^2, 1);
B = zeros((lmax+1)^2, 1);

for l = 0 : lmax
    idx = (l^2+1) : (l^2+2*l+1);
    [t_temp, A_temp, B_temp] = local_solver(l, z, R, pot_radial, N);    
    t_mat(idx, 1) = t_temp;
    A(idx) = A_temp;
    B(idx) = B_temp;    
end    
for j = 2 : numAtom
    t_mat(:, j) = t_mat(:, 1);
end

% Change some of the central local potentials manually
% remove some sites randomly
if (1)
site_idx = (numAtom+1)/2;  % the central site is of interest
rand_idx = [1 2 5 -2 -3];
vacancy_idx = site_idx - rand_idx;
for j = 1 : size(rand_idx, 2)
    t_mat(:, vacancy_idx(j)) = zeros((lmax+1)^2, 1);
end
end


T = diag(reshape(t_mat, [ ], 1));
    
% The interpolation is used when z is on the real axis
if real(z) >0 && abs(imag(z))<=1e-2
%     z1 = z + 1i;
    z1 = z + 1i*0.1;
    [G_ref1, T_ref1] = structural_green(z1, lmax, numAtom, R, offset, pot_radial_ref, N);    
%     z2 = z - 1i;
    z2 = z - 1i*0.1;
    
    [G_ref2, T_ref2] = structural_green(z2, lmax, numAtom, R, offset, pot_radial_ref, N);
    G_ref = (G_ref1+G_ref2)/2;
    T_ref = (T_ref1+T_ref2)/2;
else
    [G_ref, T_ref] = structural_green(z, lmax, numAtom, R, offset, pot_radial_ref, N);
end

% Check the decay rate of G_ref (all entries)
if (0)
x = zeros(numAtom, numAtom);
y = zeros(numAtom, numAtom);
for j = 1 : numAtom
    for n = 1 : numAtom
        x(j, n) = abs(j-n);
        y(j, n) = norm(G_ref((lmax+1)^2*(j-1)+1:(lmax+1)^2*j, ...
            (lmax+1)^2*(n-1)+1:(lmax+1)^2*n), 'fro');
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
fprintf('The 1-norm of contraction matrix is %f.\n', Q_norm1)
end


end