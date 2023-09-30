% Solve (-d^2/dx^2 + V(x)) u(x) = Eu(x) within MST 
% Check the convergence of LDoS w.r.t. the site truncation in 1D atom chain
% The M-matrix only needs to be assembled once, and the truncation can be
% directly performed on the total (large) M-matrix.
%
% 15/Sep/2023

%% Problem setting
clear
numAtom_exact =101;  
tr_largest = 40;  

CentralAtom = 0.0;   % x-coordinate of the central atom
AtomDist = 1.0;

FirstAtom = CentralAtom - AtomDist*(numAtom_exact-1)/2;

% atom positions
offset = (FirstAtom : AtomDist : FirstAtom+AtomDist*(numAtom_exact-1));  
R = AtomDist/2;   

% identical soft Coulomb potentials
if (1)
load rand_position.mat a
pot_radial = @(x) -1 ./ sqrt(1+x.^2) .* (x<=R);
pot = @(x) -1 ./ sqrt(1+(x-0.5).^2 ) .* (abs(x-0.5)<=R);   

R0 = 0.2;  
pot_radial_ref = 20;    
end


%% Reference information
% Solve the periodic eigenvalue problem and plot the band structure of reference
if (1)
LsCell = AtomDist;    % simulation cell is [0, LsCell]

% Note that the constant potential in the whole cell is meanless
% constant potential within the artificial MT
pot_ref = @(x) pot_radial_ref * (abs(x-LsCell/2)<=R0);   

NsCell_ref = LsCell * 10;    % discretization number 
numOrbital_ref = 2; 
eig_ref = eig_solver_1d(LsCell, NsCell_ref, numOrbital_ref, pot_ref, 0);
Eb_ref = min(eig_ref(1, :));           % lower bound of energy band
fprintf('The lower bound of band for reference system is %f\n', Eb_ref)
end


%% Generate the contour enclosing the desired spectrum
% Solve an inaccurate eigenvalue problem with pbc using dual planewaves
% to roughly determine the energy range, 
% then draw a coutour enclosing the desired spectrum

LsCell = AtomDist;    % simulation cell is [0, LsCell]
NsCell = LsCell * 10;    % discretization number
numOrbital = 2; 
eig = eig_solver_1d(LsCell, NsCell, numOrbital, pot, 0);
Eb = min(eig(1, :));           % lower bound of energy band
HOMO = max(eig(1, :));   % highest occupied molecular orbital energy
LUMO = min(eig(2, :));     % lowest unoccupied molecular orbital energy
gap = LUMO - HOMO;
% EF = Eb + (Eb_ref - Eb)/2;
% EF = 2;
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

% Check the energy sampling on the complex plain
if (0)
   figure 
   plot(energy_point(:, 1), energy_point(:, 2), 'bo', 'markersize', 10)
   hold on
   plot([Eb HOMO LUMO], [0 0 0], 'r*', 'markersize', 10)
   plot(Eb_ref, 0, 'c+', 'markersize', 10)
end


%% Estimate the local density of state 
site_idx = (numAtom_exact+1)/2;  % the central site is of interest
D_tr = zeros(tr_largest, 1);
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

%% Test the convergence rate of X for a fixed energy point
if (0)
err_X = zeros(tr_largest, 1);
z = energy_point(num_Matsubara, 1) + 1i * energy_point(num_Matsubara, 2);
% z = EF + 0.1*1i;
fprintf('The fixed energy point is %f+%f i\n', real(z), imag(z))
[M, A, B, T_diff] = M_basis_chain_1d(z, numAtom_exact, R, offset, pot_radial, ...
                                        R0, pot_radial_ref);
X = M^(-1);
Xn = X(2*site_idx-1 : 2*site_idx, 2*site_idx-1 : 2*site_idx);   
% submatrix corresponding to the site of interest
for dist = 1 : tr_largest
   numAtom_tr = 1 + 2*dist;
   Firstsite_idx = (numAtom_exact - numAtom_tr)/2+1;  
   Lastsite_idx = Firstsite_idx+numAtom_tr-1;
   M_tr = M(2*Firstsite_idx-1 : 2*Lastsite_idx, 2*Firstsite_idx-1 : 2*Lastsite_idx);
   X_tr = M_tr^(-1);
  
   site_idx_tr = dist+1;  % index of the central site
   Xn_tr = X_tr(2*site_idx_tr-1 : 2*site_idx_tr, 2*site_idx_tr-1 : 2*site_idx_tr);     
   err_X(dist) = norm(Xn-Xn_tr, 'fro');
end
figure
semilogy(1 : tr_largest, err_X, 'b-o',  'linewidth', 2, 'markersize', 15)
xlabel('$R$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$|X^R_{nn}(z)-X_{nn}(z)|$', 'interpreter', 'latex', 'fontsize', 20)
end


%% Note that the integration (or residue) weight are different for energy points
if (1)
D_exact = 0.0;
for j = 1 : num_energy
    z = energy_point(j, 1) + 1i * energy_point(j, 2);
    [M, A, B, T_diff] = M_alloy(a, z, numAtom_exact, R, offset, pot_radial, ...
                                        R0, pot_radial_ref);
    X = M^(-1);
    Xn = X(2*site_idx-1 : 2*site_idx, 2*site_idx-1 : 2*site_idx);   
    % submatrix corresponding to the site of interest
    Tn_diff = T_diff(2*site_idx-1 : 2*site_idx, 2*site_idx-1 : 2*site_idx); 
    Gn = -Tn_diff^(-1) + Tn_diff^(-1)*Xn*Tn_diff^(-1);
    D_exact = D_exact - weight(j)*sum(B) + weight(j)*sum(A.*diag(Gn));
    for dist = 1 : tr_largest
        numAtom_tr = 1 + 2*dist;
        Firstsite_idx = (numAtom_exact - numAtom_tr)/2+1;  
        Lastsite_idx = Firstsite_idx+numAtom_tr-1;
        M_tr = M(2*Firstsite_idx-1 : 2*Lastsite_idx, 2*Firstsite_idx-1 : 2*Lastsite_idx);
        X_tr = M_tr^(-1);
        site_idx_tr = dist+1;  % index of the central site
        Xn_tr = X_tr(2*site_idx_tr-1 : 2*site_idx_tr, 2*site_idx_tr-1 : 2*site_idx_tr);   
        Gn_tr = -Tn_diff^(-1) + Tn_diff^(-1)*Xn_tr*Tn_diff^(-1);
        D_tr(dist) = D_tr(dist) - weight(j)*sum(B) + weight(j)*sum(A.*diag(Gn_tr)); 
    end
end
D_exact = -D_exact/(2*pi*1i);
D_tr = -D_tr/(2*pi*1i);
err = abs(D_tr - D_exact);

figure
semilogy(1:40, err(1:40), '-o', 'color', [0.3010 0.7450 0.9330], ...
                    'linewidth', 2, 'markersize', 8, 'markerfacecolor', [0.3010 0.7450 0.9330])
xlabel('$R$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$|D_n^R(f)-D_n(f)|$', 'interpreter', 'latex', 'fontsize', 20)
end