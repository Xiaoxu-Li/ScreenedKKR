function eig = eig_solver_aniso(LsCell, NsCell, numOrbital, pot, Isplot)
% Solve the periodic eigenvalue problem (-\nabla + V(r)) u(r) = Eu(r)
% and plot the band structure
% Only work for the quasi-1D systems, where the k-points sampling is along one
% direction (x-direction) and choose the Gamma point for the other two directions

LsCellx = LsCell(1);
LsCelly = LsCell(2);
LsCellz = LsCell(3);
NsCellx = NsCell(1);
NsCelly = NsCell(2);
NsCellz = NsCell(3);

hs = LsCellx / NsCellx;     %   discretization size
NsCell3D = prod(NsCell);  %   number of grid points
numKpt = 20;  % number of k-points

kptGrid = -pi /LsCellx + 2*pi / LsCellx * (0:numKpt-1)'/numKpt;  
[XkptGrid3D, YkptGrid3D, ZkptGrid3D] = ndgrid(kptGrid, 0, 0);
kptGrid3D = [XkptGrid3D(:), YkptGrid3D(:), ZkptGrid3D(:)];  % k-points in 3D

% uniform grid in the real space
rptGridx = (0:NsCellx-1)' * hs;
rptGridy = (0:NsCelly-1)' * hs;
rptGridz = (0:NsCellz-1)' * hs;
[XrptGrid3D, YrptGrid3D, ZrptGrid3D] = ndgrid(rptGridx, rptGridy, rptGridz);
rptGrid3D = [XrptGrid3D(:), YrptGrid3D(:), ZrptGrid3D(:)];     

% potential field
VrptGrid3D = pot(rptGrid3D);

% uniform grid in the reciprocal space
GptGridx = (2*pi*[0:NsCellx/2, -NsCellx/2+1:-1]/LsCellx)';
GptGridy = (2*pi*[0:NsCelly/2, -NsCelly/2+1:-1]/LsCelly)';
GptGridz = (2*pi*[0:NsCellz/2, -NsCellz/2+1:-1]/LsCellz)';
[XGptGrid,YGptGrid,ZGptGrid] = ndgrid(GptGridx, GptGridy, GptGridz);
GptGrid3D = [XGptGrid(:) YGptGrid(:) ZGptGrid(:)];

eig = zeros(numOrbital, numKpt);
for k = 1 : numKpt
% kinetic term 
% (Note that there is no 1/2 before nabla !!!)
LapMult = sum((GptGrid3D+kptGrid3D(k, :)).^2, 2);
% LapMult = 0.5 * sum((GptGrid3D+kptGrid3D(k, :)).^2, 2);
LapMult = reshape(LapMult, NsCell);

% convolution via the reciprocal space
lapFunc = @(x) reshape( ...
  ifftn(LapMult .* fftn(reshape(x, NsCell))), [NsCell3D,1]);

% potential field term
VFunc = @(x) VrptGrid3D.*x;

% Fock operator in the real-space representation. 
H = @(x) lapFunc(x) + VFunc(x);

% Matrix form of the Fock operator
HMat = zeros(NsCell3D, NsCell3D);
IMat = eye(NsCell3D);
for j = 1 : NsCell3D
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
    plot(kptGrid, E_temp', 'b-o', 'linewidth', 2, 'markersize', 10)
%     plot(kptGrid, E_temp', 'b-', 'linewidth', 2, 'markersize', 12)
end
end

end
