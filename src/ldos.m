function D = ldos(lmax, numAtom)
% This script solves (-\nabla + V(r)) u(r) = Eu(r) within MST 
% and generates the local density of state (LDoS) with a test function.
% LDoS = \int_{C} f(z) \int_{\Omega_n} G(r,r;z) dr dz
%
% 28/JAN/2023


%% Problem setting
% The quasi-1D atom chain is centered at [0, 0.5, 0.5], which stretched
% along x-axis with (numAtom) atoms

CentralAtom = 0;     % x-coordinate of the central atom
AtomDist = 1.0;
FirstAtom = CentralAtom - AtomDist*(numAtom-1)/2;

offset = 0.5 * ones(numAtom, 3); % atom positions
offset(:, 1) = (FirstAtom : AtomDist : FirstAtom+AtomDist*(numAtom-1))';  
R = 0.4;   % radius of MT potential sphere

%  identical coulomb potentials
if (0) 
pot_radial = @(x) -10 ./ abs(x) + 10/R;  
pot = @(x) generate_coul(x, R, 0.5*ones(1,3), pot_radial);  % potential in each cell
end

% identical Gaussian potentials
if (1)
offset_cell = 0.5*ones(1, 3);
pot_radial = @(x) 100 * (exp(-x.^2 / 1^2) - exp(-R^2/1^2)) .* (x<=R);
pot = @(x) 100 * (exp(- sum((x-offset_cell).^2, 2) / 1^2) ...
    - exp(-R^2/1^2)) .* (sum((x-offset_cell).^2, 2)<=R^2);   % potential in each cell
end


%% Generate the contour enclosing the desired spectrum
% Solve an inaccurate eigenvalue problem with pbc using dual planewaves
% to roughly determine the energy range, 
% then draw a coutour enclosing the desired spectrum

LsCellx = AtomDist;    % simulation cell is [0, LsCellx]*[0, LsCelly]*[0, LsCellz]
LsCelly = AtomDist;
LsCellz = AtomDist;
LsCell = [LsCellx LsCelly LsCellz];
NsCell = LsCell .* [10 10 10];    % discretization number along each direction
numOrbital = 2; 
% eig = eig_solver(LsCell, NsCell, numOrbital, pot);
eig = eig_solver_aniso(LsCell, NsCell, numOrbital, pot, 0);
Eb = min(eig(1, :));           % lower bound of energy band
HOMO = max(eig(1, :));   % highest occupied molecular orbital energy
LUMO = min(eig(2, :));     % lowest unoccupied molecular orbital energy
gap = LUMO - HOMO;
EF = HOMO + 0.2*gap;
fprintf('The lower bound of band is %f\n', Eb)
fprintf('The highest occupied molecular orbital energy is %f\n', HOMO)
fprintf('The lowest unoccupied molecular orbital energy is %f and gap is %f\n', ...
    LUMO, gap)

 El = max(eig(1) - gap/2, 0.1);
Er = eig(1) + gap/2;
Gamma = 0.2 * gap;
num_energy_real = 10;  % number of energy points along real axis
num_energy_imag = 3;  % number of energy points along imaginary axis
num_energy = 2*num_energy_real+2*num_energy_imag;
energy_point = zeros(num_energy, 2);
% first column for real part and second column for imaginary part

% construct a rectangular contour [El, Er] * [-Gamma, Gamma]
length_contour = 4*Gamma + 2*(Er-El);
% left side
for j = 1 : num_energy_imag
    energy_point(j, 2) = -Gamma + 2*Gamma / (num_energy_imag+1) * j;
    energy_point(j, 1) = El;
end
% right side
energy_point(num_energy_imag+1 : 2*num_energy_imag, 1) = Er;
energy_point(num_energy_imag+1 : 2*num_energy_imag, 2) ...
    = energy_point(1 : num_energy_imag, 2);
% top
for j = 1 : num_energy_real
    energy_point(j+2*num_energy_imag, 1) = El + (Er-El) / (num_energy_real-1) * (j-1);
    energy_point(j+2*num_energy_imag, 2) = Gamma;
end
% bottom
energy_point(2*num_energy_imag+num_energy_real+1 : end, 1) ...
    = energy_point(2*num_energy_imag+1 : 2*num_energy_imag+num_energy_real, 1) ;
energy_point(2*num_energy_imag+num_energy_real+1 : end, 2) = -Gamma;


%% Estimate the local density of state 
site_idx = (numAtom+1)/2;  % the central site is of interest
D = 0.0;
for j = 1 : num_energy
    z = energy_point(j, 1) + 1i * energy_point(j, 2);
    [M, A, B, T_diff] = M_basis_chain(z, lmax, numAtom, R, offset, pot_radial);
    X = M^(-1);
    Xn = X((lmax+1)^2*(site_idx-1)+1 : (lmax+1)^2*site_idx, ...
        (lmax+1)^2*(site_idx-1)+1 : (lmax+1)^2*site_idx);   
    Tn_diff = T_diff((lmax+1)^2*(site_idx-1)+1 : (lmax+1)^2*site_idx, ...
        (lmax+1)^2*(site_idx-1)+1 : (lmax+1)^2*site_idx); 
    Gn = -Tn_diff^(-1) + Tn_diff^(-1)*Xn*Tn_diff^(-1);
    % submatrix corresponding to the site of interest
    D = D - fermi_dirac(1,EF,z)*sum(B) + fermi_dirac(1,EF,z)*sum(A.*diag(Gn));    
end
D = length_contour/num_energy * D;

end
