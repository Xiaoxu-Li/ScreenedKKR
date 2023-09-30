 function M = M_chain(z, lmax, numAtom, R, offset, pot_radial)
% Consider (-\nabla + V(r)) u(r) = Eu(r) within MST and
% Generate the M-matrix for an atom chain
% 
% 02/FEB/2023

N = 1e2;
t_mat = zeros((lmax+1)^2, numAtom);     % t-matrix
for l = 0 : lmax
    idx = (l^2+1) : (l^2+2*l+1);
    t_mat(idx, 1) = t_matrix_fd(l, z, R, pot_radial, N);    
end    
for j = 2 : numAtom
    t_mat(:, j) = t_mat(:, 1);
end
T = diag(reshape(t_mat, [ ], 1));

% pot_radial_ref = 50;    % constant potential within MT
pot_radial_ref = 120; 

% The interpolation is used when z is on the real axis
% Note that this operation has a drawback. 
% For the energy in the spectrum of reference, G_ref has an incorrect exponential
% off-diagonal decay, which leads to the wrong decay property of X.
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
M = T_diff^(-1) - G_ref;   % M-matrix used to be inverted

end