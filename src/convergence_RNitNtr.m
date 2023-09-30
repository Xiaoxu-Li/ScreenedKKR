function [err_it, err_tr]= convergence_RNitNtr(lmax)
% Solve (-\nabla + V(r)) u(r) = Eu(r) within MST 
% Check the convergence of LDoS with respect to Nit and Ntr
% in a quasi-1D atom chain.
%
% 09/Aug/2023

%% Problem setting
tic
% addpath(genpath('/Users/xxl/Desktop/Matlab-root-tracing'));

% lmax = 0;
numAtom_exact = 51;

CentralAtom = 0.0;   % x-coordinate of the central atom
AtomDist = 1.0;
FirstAtom = CentralAtom - AtomDist*(numAtom_exact-1)/2;

offset = 0.5 * ones(numAtom_exact, 3); % atom positions
offset(:, 1) = (FirstAtom : AtomDist : FirstAtom+AtomDist*(numAtom_exact-1))';  

% identical Gaussian potentials
if (1)
R = 0.3;   % radius of MT potential sphere

center_offset = 0.5 * ones(1, 3);
pot_radial = @(x) 10 * (exp(- x.^2 / 0.1^2) - exp(-R^2/0.1^2)) .* (x<=R);
pot = @(x) 10 * exp(- sum((x-center_offset).^2, 2) /0.1^2) ...
            .* (sum((x-center_offset).^2, 2)<=R^2);    % potential in each cell
pot_radial_ref = 10;    % constant potential within MT
end


%% Reference information
% Solve the periodic eigenvalue problem and plot the band structure of reference
if (1)
LsCellx = AtomDist;    % simulation cell is [0, LsCellx]*[0, LsCelly]*[0, LsCellz]
LsCelly = AtomDist;
LsCellz = AtomDist;
LsCell = [LsCellx LsCelly LsCellz];
NsCell = LsCell .* [10 10 10];    % discretization number along each direction

pot_ref = @(x) pot_radial_ref * (sum((x-center_offset).^2, 2)<=R^2);   

numOrbital_ref = 2; 
eig_ref = eig_solver_aniso(LsCell, NsCell, numOrbital_ref, pot_ref, 0);
Eb_ref = min(eig_ref(1, :));           % lower bound of energy band
fprintf('The lower bound of band for reference system is %f\n', Eb_ref)
end


%% Generate the contour enclosing the desired spectrum
% Solve an inaccurate eigenvalue problem with pbc using dual planewaves
% to roughly determine the energy range, 
% then draw a coutour enclosing the desired spectrum

numOrbital = 2; 
eig = eig_solver_aniso(LsCell, NsCell, numOrbital, pot, 0);
Eb = min(eig(1, :));           % lower bound of energy band
HOMO = max(eig(1, :));   % highest occupied molecular orbital energy
LUMO = min(eig(2, :));     % lowest unoccupied molecular orbital energy
gap = LUMO - HOMO;
% EF = Eb + (Eb_ref - Eb)/2;
EF = 0.5;

fprintf('The lower bound of band is %f\n', Eb)
fprintf('The highest occupied molecular orbital energy is %f\n', HOMO)
fprintf('The lowest unoccupied molecular orbital energy is %f and gap is %f\n', ...
    LUMO, gap)
fprintf('The fixed Fermi level is %f\n', EF)

El = min(-0.5, Eb-1);
Er = HOMO + gap/2;
beta = 2;      % beta=1/(k*T)
num_Matsubara = 5;    
% half number (due to symmetry) of Matsabara points encloed in the contour
Gamma = 2*num_Matsubara*pi/beta;   
% Gamma should be chosen to avoid the poles of Fermi Dirac function
% (Matsubara points) on the imaginary axis, i.e. (2*j+1)*pi*k*T=(2*j+1)*pi/beta

num_energy_real = 10;  % number of energy points along real axis
num_energy_imag = 5;  % number of energy points along imaginary axis
num_energy = 2*num_energy_real+num_energy_imag + 2*num_Matsubara;
energy_point = zeros(num_energy, 2);
% first column for real part and second column for imaginary part

% Construct a rectangular contour [El, Er] * [-Gamma, Gamma] without right
% side, the points on the right side can be omited due to fast decay
% of fermi-dirac function with low-temperature
% Note that the contour is clockwise!

% Matsabara points
for j = 1 : 2*num_Matsubara
    energy_point(j, 2) = -(2*num_Matsubara-1)*pi/beta + 2*(j-1)*pi/beta;
    energy_point(j, 1) = EF;
end
% left side
for j = 1 : num_energy_imag
    energy_point(2*num_Matsubara+j, 2) = ...
        -Gamma + (j-1)*2*Gamma / (num_energy_imag-1);
    energy_point(2*num_Matsubara+j, 1) = El;
end
% top
for j = 1 : num_energy_real
    energy_point(2*num_Matsubara+num_energy_imag+j, 1) ...
        = El + j*(Er-El) / num_energy_real;
    energy_point(2*num_Matsubara+num_energy_imag+j, 2) = Gamma;
end
% bottom
energy_point(2*num_Matsubara+num_energy_imag+num_energy_real+1 : end, 1) ...
    = energy_point(2*num_Matsubara+num_energy_imag+1 : ...
                            2*num_Matsubara+num_energy_imag+num_energy_real, 1) ;
energy_point(2*num_Matsubara+num_energy_imag+num_energy_real+1 : end, 2) = -Gamma;


%% Estimate the local density of state 
site_idx = (numAtom_exact+1)/2;  % the central site is of interest
D_exact = 0.0;
% NtrMax = (numAtom_exact-1)/2;   % truncation for G_ref
NtrMax = 20;
NitMax = 20;    % the maximum iteration number
initial = zeros((lmax+1)^2*numAtom_exact, (lmax+1)^2);
initial((lmax+1)^2*(site_idx-1)+1 : (lmax+1)^2*site_idx, :) = eye((lmax+1)^2);
% initial guess for n-th column of X
weight = zeros(num_energy, 1);

% Calculate the weight matrix for energy numerical integration
weight(1: 2*num_Matsubara) = -2*pi*1i/beta;   % residue
% left
for j = 1 : num_energy_imag
    z = energy_point(2*num_Matsubara+j, 1) ...
                + 1i * energy_point(2*num_Matsubara+j, 2);
    weight(2*num_Matsubara+j) = (2*Gamma*1i)/num_energy_imag...
                    *fermi_dirac(beta,EF,z);
end
% top
for j = 1 : num_energy_real
    z = energy_point(2*num_Matsubara+num_energy_imag+j, 1) ...
        + 1i * energy_point(2*num_Matsubara+num_energy_imag+j, 2);
     weight(2*num_Matsubara+num_energy_imag+j) ...
         = (Er-El)/num_energy_real*fermi_dirac(beta,EF,z);
end
% bottom
for j = 1 : num_energy_real
   z = energy_point(2*num_Matsubara+num_energy_imag+num_energy_real+j, 1) ...
        + 1i * energy_point(2*num_Matsubara+num_energy_imag+num_energy_real+j, 2);
     weight(2*num_Matsubara+num_energy_imag+num_energy_real+j) ...
         = - (Er-El)/num_energy_real*fermi_dirac(beta,EF,z);
end

% Test the convergence rate for LDOS
if (1)
D_it = zeros(NitMax, 1);
D_tr = zeros(NtrMax, 1);
for j = 1 : num_energy
    z = energy_point(j, 1) + 1i * energy_point(j, 2);
    [G_ref, A, B, T_diff] = Gref_basis_chain(z, lmax, numAtom_exact, R, offset, ...
                                            pot_radial, pot_radial_ref);
    M = T_diff^(-1) - G_ref;   
    X = M^(-1);
    Xn = X((lmax+1)^2*(site_idx-1)+1 : (lmax+1)^2*site_idx, ...
        (lmax+1)^2*(site_idx-1)+1 : (lmax+1)^2*site_idx);   
    % submatrix corresponding to the site of interest
    Tn_diff = T_diff((lmax+1)^2*(site_idx-1)+1 : (lmax+1)^2*site_idx, ...
        (lmax+1)^2*(site_idx-1)+1 : (lmax+1)^2*site_idx); 
    Gn = -Tn_diff^(-1) + Tn_diff^(-1)*Xn*Tn_diff^(-1);
    D_exact = D_exact - weight(j)*sum(B) + weight(j)*sum(A.*diag(Gn)); 
    
    % fix Ntr, check the convergence w.r.t. Nit
    for Nit = 1 : NitMax
        X_tr = tfqmr_3d(NtrMax, Nit, initial, lmax, site_idx, G_ref, T_diff);
        Xn_tr = X_tr((lmax+1)^2*(site_idx-1)+1 : (lmax+1)^2*site_idx, :);
        Gn_tr = -Tn_diff^(-1) + Tn_diff^(-1)*Xn_tr*Tn_diff^(-1);
        D_it(Nit) = D_it(Nit) - weight(j)*sum(B) + weight(j)*sum(A.*diag(Gn_tr));
    end

    % fix Nit, check the convergence w.r.t. Ntr
    for Ntr = 1 : NtrMax
        X_tr = tfqmr_3d(Ntr, NitMax, initial, lmax, site_idx, G_ref, T_diff);
        Xn_tr = X_tr((lmax+1)^2*(site_idx-1)+1 : (lmax+1)^2*site_idx, :);
        Gn_tr = -Tn_diff^(-1) + Tn_diff^(-1)*Xn_tr*Tn_diff^(-1);
        D_tr(Ntr) = D_tr(Ntr) - weight(j)*sum(B) + weight(j)*sum(A.*diag(Gn_tr));
    end
end
D_exact = -D_exact/(2*pi*1i);
D_tr = -D_tr/(2*pi*1i);
D_it = -D_it/(2*pi*1i);
err_tr = abs(D_tr - D_exact);
err_it = abs(D_it - D_exact);

% figure
% semilogy(1:NitMax, err_it, 'b-d',  'linewidth', 2, 'markersize', 20)
% xlabel('$N_{it}$', 'interpreter', 'latex', 'fontsize', 20)
% ylabel('$|D_n^{RN_{it}N_{tr}}(f)-D_n(f)|$', 'interpreter', 'latex', 'fontsize', 20)
% title('Fix $R$ and $N_{tr}$', 'interpreter', 'latex', 'fontsize', 20)
% 
% figure
% semilogy(1:NtrMax, err_tr, 'b-d',  'linewidth', 2, 'markersize', 20)
% xlabel('$N_{tr}$', 'interpreter', 'latex', 'fontsize', 20)
% ylabel('$|D_n^{RN_{it}N_{tr}}(f)-D_n(f)|$', 'interpreter', 'latex', 'fontsize', 20)
% title('Fix $R$ and $N_{it}$', 'interpreter', 'latex', 'fontsize', 20)

toc

end