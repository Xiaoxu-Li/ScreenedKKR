function [G_ref, T_ref] = structural_green(z, lmax, numAtom, R, offset, pot_radial_ref, N)
% Compute the structural Green's function and T-matrix of the reference
% for a given energy parameter z
% 
% 27/JAN/2023

%% Solve the periodic eigenvalue problem 
% and plot the band structure of reference
if (0)
offset_cell = 0.5*ones(1, 3); 
% constant potential within MT
pot = @(x) pot_radial_ref * (sum((x-offset_cell).^2, 2)<=R^2);   

AtomDist = 1.0;
LsCellx = AtomDist;    % simulation cell is [0, LsCellx]*[0, LsCelly]*[0, LsCellz]
LsCelly = AtomDist;
LsCellz = AtomDist;
LsCell = [LsCellx LsCelly LsCellz];
NsCell = LsCell .* [4 4 4];    % discretization number along each direction
numOrbital = 2; 
eig = eig_solver_aniso(LsCell, NsCell, numOrbital, pot, 0);
Eb = min(eig(1, :));           % lower bound of energy band
fprintf('The lower bound of band for reference system is %f\n', Eb)
end

%% structure constant
g = zeros((lmax+1)^2, numAtom, (lmax+1)^2, numAtom);   
for l1 = 0 : lmax
    for m1 = -l1 : l1
        idx1 = l1^2+l1+m1+1;            
        for l2 = 0 : lmax
            for m2 = -l2 : l2              
                idx2 = l2^2+l2+m2+1; 
                for n1 = 1 : numAtom
                    offset1 = offset(n1, :);
                    for n2 = 1 : numAtom
                        if n1==n2
                            continue
                        else
                            offset2 = offset(n2, :);
                            g(idx2, n2, idx1, n1) = structure_const(l1,l2,m1,m2, z, 10, offset2-offset1);
                        % using larger cut-off when computing structure constants
                        end
                    end
                end
            end
        end
    end
end
G0 = reshape(g, numAtom*(lmax+1)^2, numAtom*(lmax+1)^2);    

%% Structural Green's function of reference
shift_z = z - pot_radial_ref;
h = R/N; 
xgrid = (1:N)'*h;

% Guarantee the imaginary part of sqrt(z) is positive
if imag(sqrt(z))>=0
    sqrtz = sqrt(z); 
else
    sqrtz = -sqrt(z);
end

if imag(sqrt(shift_z))>=0
    sqrtzShift = sqrt(shift_z); 
else
    sqrtzShift = -sqrt(shift_z);
end

t_mat_ref = zeros((lmax+1)^2, numAtom);     
for l = 0 : lmax
    idx = (l^2+1) : (l^2+2*l+1);
    tl = pot_radial_ref * sum( spherical_bessel(l, sqrtz*xgrid) .* xgrid .^2 ...
        .* spherical_bessel(l, sqrtzShift*xgrid)) * h;
    tl = -1i * sqrtz * tl;
    normalization = spherical_bessel(l, sqrtz*R) / (1-tl*spherical_hankel(l, sqrtz*R));
    tl = tl * normalization;
    t_mat_ref(idx, 1) = tl; 
end    
for j = 2 : numAtom
    t_mat_ref(:, j) = t_mat_ref(:, 1);
end
T_ref = diag(reshape(t_mat_ref, [ ], 1));

G_ref = (eye((lmax+1)^2*numAtom) - G0*T_ref)^(-1) * G0;

% Check the decay rate of G_ref (all entries)
if (0)
x = zeros(numAtom, numAtom);
y = zeros(numAtom, numAtom);
for j = 1 : numAtom
    for n = 1 : numAtom
        x(j, n) = abs(j-n);
        y(j, n) = norm(G_ref((lmax+1)^2*(j-1)+1:(lmax+1)^2*j, ...
            (lmax+1)^2*(n-1)+1:(lmax+1)^2*n), 'fro');
    end
end
x = reshape(x, [ ], 1);
y = reshape(y, [ ], 1);
figure
semilogy(x, y, 'bo', 'markersize', 10)
end