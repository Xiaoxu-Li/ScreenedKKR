function X = tfqmr_1d(Ntr, Nit, initial, site_idx, G_ref, T_diff)
% Using transpose-free quasi-minimal residual scheme 
% to compute X instead of directly invertion of M
% for 1D system (i.e. L=2)
%
% Output: X is a column vector corresponding to n-th site
% 08/Aug/2023

% truncated structural Green's function of reference
numAtom = size(G_ref, 1)/2;
blk = ones(2);
B = zeros(numAtom, numAtom);
for j = 1 : Ntr
   B = B + diag(ones(numAtom-j,1), j);
end
B = B + B';
B = B + eye(numAtom);
indicator = kron(B, blk);
G_ref_tr = G_ref .* indicator;

% solve AX = b
Q = T_diff * G_ref_tr;
A = eye(2*numAtom) - Q;

b = T_diff(:, 2*(site_idx-1)+1 : 2*site_idx);
[x1,~] = tfqmr(A, b(:, 1), 1e-15, Nit, [ ], [ ], initial(:, 1));
[x2,~]= tfqmr(A, b(:, 2), 1e-15, Nit, [ ], [ ], initial(:, 2));
X = [x1 x2];

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
