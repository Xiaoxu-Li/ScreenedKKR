 function M = M_chain_1d(z, numAtom, R, offset, pot_radial)
% Consider (-d^2/dx^2 + V(x)) u(x) = Eu(x) within 1D MST and
% Generate the M-matrix for an atom chain
% 
% 29/MAR/2023

N = 1e2;
t_mat = zeros(2, numAtom);     % t-matrix
% Note that the first index is set to be l=1 and the second is l=-1.
t_mat(1,1) = t_matrix_fd_1d(1, z, R, pot_radial, N); 
t_mat(2,1) = t_matrix_fd_1d(-1, z, R, pot_radial, N); 

for j = 2 : numAtom
    t_mat(:, j) = t_mat(:, 1);
end
T = diag(reshape(t_mat, [ ], 1));

pot_radial_ref = 120;    % constant potential within the artificial MT

% The interpolation is used when z is on the real axis
% Note that this operation has a drawback. 
% For the energy in the spectrum of reference, G_ref has an incorrect exponential
% off-diagonal decay, which leads to the wrong decay property of X.
if real(z) >0 && abs(imag(z))<=1e-2
%     z1 = z + 1i;
    z1 = z + 1i*0.1;
    [G_ref1, T_ref1] = structural_green_1d(z1, numAtom, offset, pot_radial_ref, N);    
%     z2 = z - 1i;
    z2 = z - 1i*0.1;    
    [G_ref2, T_ref2] = structural_green_1d(z2, numAtom, offset, pot_radial_ref, N);
    G_ref = (G_ref1+G_ref2)/2;
    T_ref = (T_ref1+T_ref2)/2;
else
    [G_ref, T_ref] = structural_green_1d(z, numAtom, offset, pot_radial_ref, N);
end

T_diff = T - T_ref;
% cond(T_diff)
M = T_diff^(-1) - G_ref;   % M-matrix used to be inverted

end