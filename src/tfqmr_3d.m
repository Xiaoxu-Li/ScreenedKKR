function X = tfqmr_3d(Ntr, Nit, initial, lmax, site_idx, G_ref, T_diff)
% Using transpose-free quasi-minimal residual scheme 
% to compute X instead of directly invertion of M
%
% Output: X is a column vector corresponding to n-th site
%
% 09/Aug/2023

% truncated structural Green's function of reference
L = (lmax+1)^2;
numAtom = size(G_ref, 1)/L;
blk = ones(L);
A = zeros(numAtom, numAtom);
for j = 1 : Ntr
   A = A + diag(ones(numAtom-j,1), j);
end
A = A + A';
A = A + eye(numAtom);
indicator = kron(A, blk);
G_ref_tr = G_ref .* indicator;

% solve AX = b
Q = T_diff * G_ref_tr;
A = eye(L*numAtom) - Q;
X = zeros(L*numAtom, L);
b = T_diff(:, L*(site_idx-1)+1 : L*site_idx);
for j = 1 : L
    [x, ~] = tfqmr(A, b(:, j), 1e-15, Nit, [ ], [ ], initial(:, j));
    X(:, j) = x;
end


% Test the convergence of iterative scheme
if (0)
NitMax = 50;
X_exact = (eye(size(Q, 1))-Q)^(-1) * Tn_diff;
err = zeros(NitMax, 1);

for Nit = 1 : NitMax
    X = initial;
    for k = 1 : Nit
        X = Tn_diff + Q*X;
    end
    err(Nit) = abs(X(site_idx)-X_exact(site_idx, site_idx));
end
figure
semilogy(1:NitMax, err, 'b-o',  'linewidth', 2, 'markersize', 12)
end
