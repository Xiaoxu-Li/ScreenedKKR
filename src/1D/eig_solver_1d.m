function eig = eig_solver_1d(LsCell, NsCell, numOrbital, pot, Isplot)
% Solve 1D periodic eigenvalue problem (-d^2/dx^2+ V(x)) u(x) = Eu(x) 
% and plot the band structure
%
% 23/MAR/2023

hs = LsCell / NsCell;     %   discretization size
numKpt = 20;  % number of k-points

kptGrid = -pi /LsCell + 2*pi / LsCell * (0:numKpt-1)'/numKpt;  % k-points
rptGrid = (0:NsCell-1)' * hs;  % uniform grid in the real space
VrptGrid = pot(rptGrid);  % potential field
GptGrid = (2*pi*[0:NsCell/2, -NsCell/2+1:-1]/LsCell)';
% uniform grid in the reciprocal space

eig = zeros(numOrbital, numKpt);
for k = 1 : numKpt
% kinetic term 
% (Note that there is no 1/2 before nabla !!!)
LapMult = (GptGrid+kptGrid(k)).^2;

% convolution via the reciprocal space
lapFunc = @(x) ifftn(LapMult .* fftn(x));

% potential field term
VFunc = @(x) VrptGrid.*x;

% Fock operator in the real-space representation. 
H = @(x) lapFunc(x) + VFunc(x);

% Matrix form of the Fock operator
HMat = zeros(NsCell, NsCell);
IMat = eye(NsCell);
for j = 1 : NsCell
  HMat(:, j) = H(IMat(:, j));
end
HMat = 0.5 * (HMat + HMat');   

[~, E] = eigs(HMat, numOrbital, 'sr');
eig(:, k) = real(diag(E));
end

% plot the band structure
if (Isplot)
figure
hold on
for j = 1 : numOrbital
    E_temp = eig(j, :);
    plot(kptGrid, E_temp', 'b-', 'linewidth', 2, 'markersize', 12)
end
end

end
