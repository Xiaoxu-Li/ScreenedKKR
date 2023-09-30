function X = mst_structural_green_1d(z)
% This script considers (-d^2/dx^2 + V(x)) u(x) = zu(x) within 1D MST 
% and generates the key quantity of structural Green's function (X-matrix)
% with a given energy parameter z.
% Refer to the manuscript for the definition of X.
%
% 29/MAR/2023

%% Set up an atom chain with several identical atoms
numAtom = 35;
AtomDist = 1.0;
FirstAtom = 0.5;
offset = (FirstAtom : AtomDist : FirstAtom+AtomDist*(numAtom-1))';  
% atom positions
R = AtomDist/2;   % no gap between two adjacent cells
center_ind = (numAtom+1)/2;

% pot_radial = @(x) -10 ./ abs(x) + 10/R; 
pot_radial = @(x) 10 * exp(- x.^2 / 1^2) .* (x<=R);
% set all the radial function to be the same


%% Solve the periodic eigenvalue problem and plot the band structure
if (1)
% potential in each cell
pot = @(x) 10 * exp(- abs(x-0.5)/ 1^2).* (abs(x-0.5)<=R);   
% pot = @(x) 100*(1./sqrt(sum((x-offset_cell).^2, 2)+1)  - 1/sqrt(R^2+1)) ...
%     .* (sum((x-offset_cell).^2, 2)<=R^2);     % smoothened Coulomb

LsCell = AtomDist;    % simulation cell is [0, LsCell]
NsCell = LsCell * 10;    % discretization number
numOrbital = 2; 
eig = eig_solver_1d(LsCell, NsCell, numOrbital, pot, 0);
Eb = min(eig(1, :));           % lower bound of energy band
HOMO = max(eig(1, :));   % highest occupied molecular orbital energy
LUMO = min(eig(2, :));     % lowest unoccupied molecular orbital energy
gap = LUMO - HOMO;
fprintf('The lower bound of band is %f\n', Eb)
fprintf('The highest occupied molecular orbital energy is %f\n', HOMO)
fprintf('The lowest unoccupied molecular orbital energy is %f and gap is %f\n', ...
    LUMO, gap)
end

%% Generate X-matrix for a given energy parameter 
M = M_chain_1d(z, numAtom, R, offset, pot_radial);
X = M^(-1);
cond(M)
eigs(M, 1, 'smallestabs')

% Check the decay rate of X (all entries)
if (1)
x = zeros(numAtom, numAtom);
y = zeros(numAtom, numAtom);
for j = 1 : numAtom
    for n = 1 : numAtom
        x(j, n) = abs(j-n);
        y(j, n) = norm(M(2*(j-1)+1:2*j, 2*(n-1)+1:2*n), 'fro');
    end
end
x = reshape(x, [ ], 1);
y = reshape(y, [ ], 1);
figure
semilogy(x, y, 'bo', 'markersize', 10, 'linewidth', 2)
end

if (1)
x = zeros(numAtom, numAtom);
y = zeros(numAtom, numAtom);
for j = 1 : numAtom
    for n = 1 : numAtom
        x(j, n) = abs(j-n);
        y(j, n) = norm(X(2*(j-1)+1:2*j, 2*(n-1)+1:2*n), 'fro');
    end
end
x = reshape(x, [ ], 1);
y = reshape(y, [ ], 1);
figure
semilogy(x, y, 'bo', 'markersize', 10, 'linewidth', 2)
end

% Check the decay rate of X (first block column)
if (0)
inv_block_norm = zeros(numAtom, 1);  % norm of first column of blocks
for n = 1 : numAtom
    inv_block_norm(n) = norm(X((n-1)*2+1 : n*2, 1:2), 'fro');
end
figure
semilogy(offset-offset(1), inv_block_norm, 'bo',  'linewidth', 2, 'markersize', 12)
end


% Calculate the truncated X for a given energy parameter and test convergence rate
if (1)
X_center = X((center_ind-1)*2+1 : center_ind*2, (center_ind-1)*2+1 : center_ind*2);

tr_largest =10;
err = zeros(tr_largest, 1);
for dist = 1 : tr_largest
    numAtom_tr = 1 + 2*dist;
    Firstind_tr = (numAtom-numAtom_tr)/2 + 1;
    Lastind_tr = Firstind_tr + numAtom_tr - 1;
    M_tr = M((Firstind_tr-1)*2+1 : Lastind_tr*2, (Firstind_tr-1)*2+1 : Lastind_tr*2);
    X_tr = M_tr^(-1);    
    center_tr_ind = dist+1;
    err(dist) = norm( X_center - X_tr((center_tr_ind-1)*2+1 : center_tr_ind*2, ...
             (center_tr_ind-1)*2+1 : center_tr_ind*2), 'fro');
end

figure
semilogy(1 : tr_largest, err, 'b-o',  'linewidth', 2, 'markersize', 20)
xlabel('$R$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$|X^R_{nn}-X_{nn}|$', 'interpreter', 'latex', 'fontsize', 20)
% title('$z=11.0$', 'interpreter', 'latex', 'fontsize', 20)
end
