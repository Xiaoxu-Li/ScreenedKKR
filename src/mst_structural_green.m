function X = mst_structural_green(z, lmax)
% This script considers (-\nabla + V(r)) u(r) = zu(r) within MST 
% and generates the key quantity of structural Green's function (X-matrix)
% with a given energy parameter z.
% Refer to the manuscript for the definition of X.
%
% 27/JAN/2023

%% Set up an atom chain with several identical atoms along x-direction
numAtom = 35;
AtomDist = 1.0;
FirstAtom = 0.5;

offset = 0.5 * ones(numAtom, 3); % atom positions
offset(:, 1) = (FirstAtom : AtomDist : FirstAtom+AtomDist*(numAtom-1))';  
R = 0.4;   % radius of MT potential sphere
center_ind = (numAtom+1)/2;

% pot_radial = @(x) -10 ./ abs(x) + 10/R; 
pot_radial = @(x) 100 * (exp(- x.^2 / 1^2) - exp(-R^2/1^2)) .* (x<=R);
% set all the radial function to be the same


%% Solve the periodic eigenvalue problem and plot the band structure
if (1)
offset_cell = 0.5*ones(1, 3); 
% potential in each cell
pot = @(x) 100 * (exp(- sum((x-offset_cell).^2, 2) / 1^2) ...
    - exp(-R^2/1^2)) .* (sum((x-offset_cell).^2, 2)<=R^2);   
% pot = @(x) 100*(1./sqrt(sum((x-offset_cell).^2, 2)+1)  - 1/sqrt(R^2+1)) ...
%     .* (sum((x-offset_cell).^2, 2)<=R^2);     % smoothened Coulomb

LsCellx = AtomDist;    % simulation cell is [0, LsCellx]*[0, LsCelly]*[0, LsCellz]
LsCelly = AtomDist;
LsCellz = AtomDist;
LsCell = [LsCellx LsCelly LsCellz];
NsCell = LsCell .* [10 10 10];    % discretization number along each direction
numOrbital = 2; 
% eig = eig_solver(LsCell, NsCell, numOrbital, pot);
eig = eig_solver_aniso(LsCell, NsCell, numOrbital, pot, 1);
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
M = M_chain(z, lmax, numAtom, R, offset, pot_radial);
X = M^(-1);

% Check the decay rate of X (all entries)
if (0)
x = zeros(numAtom, numAtom);
y = zeros(numAtom, numAtom);
for j = 1 : numAtom
    for n = 1 : numAtom
        x(j, n) = abs(j-n);
        y(j, n) = norm(X((lmax+1)^2*(j-1)+1:(lmax+1)^2*j, ...
            (lmax+1)^2*(n-1)+1:(lmax+1)^2*n), 'fro');
    end
end
x = reshape(x, [ ], 1);
y = reshape(y, [ ], 1);
figure
semilogy(x, y, 'b-o', 'markersize', 10)
end


if (0)
% Check the decay rate of X
inv_block_norm = zeros(numAtom, 1);  % norm of first column of blocks
for n = 1 : numAtom
    inv_block_norm(n) = norm(X((n-1)*(lmax+1)^2+1 : n*(lmax+1)^2, ...
        1 : (lmax+1)^2), 'fro');
end

figure
if real(z)>0 && imag(z)==0
% algebraic convergence rate for positive energies
    loglog(offset(:, 1)-offset(1, 1), inv_block_norm, 'b-o',  'linewidth', 2, 'markersize', 12)
else
% exponential convergence rate for negative/complex energies
    semilogy(offset(:, 1)-offset(1, 1), inv_block_norm, 'b-o',  'linewidth', 2, 'markersize', 12)
end
end


% Calculate the truncated X for a given energy parameter and test convergence rate
if (1)
X_center = X((center_ind-1)*(lmax+1)^2+1 : center_ind*(lmax+1)^2, ...
        (center_ind-1)*(lmax+1)^2+1 : center_ind*(lmax+1)^2);

tr_largest =10;
err = zeros(tr_largest, 1);
for dist = 1 : tr_largest
    numAtom_tr = 1 + 2*dist;
%     offset_tr = offset(center_ind-dist : center_ind+dist, :);
%     M_tr = spm_chain(z, lmax, numAtom_tr, R, offset_tr, pot_radial);    
    Firstind_tr = (numAtom-numAtom_tr)/2 + 1;
    Lastind_tr = Firstind_tr + numAtom_tr - 1;
    M_tr = M((Firstind_tr-1)*(lmax+1)^2+1 : Lastind_tr*(lmax+1)^2, ...
        (Firstind_tr-1)*(lmax+1)^2+1 : Lastind_tr*(lmax+1)^2);
    X_tr = M_tr^(-1);    
    center_tr_ind = dist+1;
    err(dist) = norm( X_center - X_tr((center_tr_ind-1)*(lmax+1)^2+1 : center_tr_ind*(lmax+1)^2, ...
             (center_tr_ind-1)*(lmax+1)^2+1 : center_tr_ind*(lmax+1)^2), 'fro');
end

figure
semilogy(1 : tr_largest, err, 'b-o',  'linewidth', 2, 'markersize', 20)
xlabel('$R$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$|X^R_{nn}-X_{nn}|$', 'interpreter', 'latex', 'fontsize', 20)
title('$z=11.0$', 'interpreter', 'latex', 'fontsize', 20)

end
