function X = Geo_expansion_1d(Ntr, P, G_ref, T_diff)
% Using geometric series to compute X instead of directly invertion of M
% for 1D system
%
% 11/July/2023

% truncated structural Green's function of reference
numAtom = size(G_ref, 1)/2;
blk = ones(2);
A = zeros(numAtom, numAtom);
for j = 1 : Ntr
   A = A + diag(ones(numAtom-j,1), j);
end
A = A + A';
A = A + eye(numAtom);
indicator = kron(A, blk);
G_ref_tr = G_ref .* indicator;

Q = G_ref_tr * T_diff;

X= eye(size(G_ref, 1)); 
for j = 1 : P-1
    X = X + Q^j;
end

% Check the assumption \|G_ref T_diff\|<1
if (0)
Q_norm1 = norm(Q, 1);
Q_norm2 = norm(Q, 2);
Q_norminf = norm(Q, inf);
fprintf('The 1, 2, inf-norm of contraction matrix is %f, %f and %f, respectively.\n', ...
    Q_norm1, Q_norm2, Q_norminf)
end
