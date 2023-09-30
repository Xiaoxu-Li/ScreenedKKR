function [G_ref, T_ref] = structural_green_1d(z, numAtom, offset, pot_radial_ref, R0, N)
% Compute the structural Green's function and T-matrix of the reference
% for a given energy parameter z in 1d system
% 
% 13/Sep/2023


%% structure constants
% Guarantee the imaginary part of sqrt(z) is positive
if imag(sqrt(z))>=0
    sqrtz = sqrt(z); 
else
    sqrtz = -sqrt(z);
end

block = eye(2);
g = zeros(2, numAtom, 2, numAtom);   
for n1 = 1 : numAtom
    for n2 = 1 : numAtom
        if n1 == n2
            continue
        else
            block(1, 2) = 1i*sign(offset(n2)-offset(n1));
            block(2, 1) = -1i*sign(offset(n2)-offset(n1));
            g(:, n1, :, n2) = exp(1i*sqrtz*abs(offset(n2)-offset(n1))) * block;
        end
    end
end
G0 = reshape(g, numAtom*2, numAtom*2);    

%% Structural Green's function of reference
shift_z = z - pot_radial_ref;
h = R0/N; 
xgrid = (1:N)'*h;

if imag(sqrt(shift_z))>=0
    sqrtzShift = sqrt(shift_z); 
else
    sqrtzShift = -sqrt(shift_z);
end

t_mat_ref = zeros(2, numAtom);     

% Note that the first index is set to be l=1 and the second is l=-1.
% t-matrix for the even component (l=1)
bessel_1d = cos(sqrtz*xgrid);
u = cos(sqrtzShift*xgrid);
tl = pot_radial_ref * sum( bessel_1d .* u) * h; 
tl = tl / (1i * sqrtz);
normalization = cos(sqrtz*R0) / (cos(sqrtzShift*R0)-tl*exp(1i*sqrtz*R0));  
t_mat_ref(1, 1) = tl * normalization;

% t-matrix for the odd component (l=-1)
bessel_1d = sin(sqrtz*xgrid);
u = sin(sqrtzShift*xgrid);
tl = pot_radial_ref * sum( bessel_1d .* u) * h; 
tl = tl / (1i * sqrtz);
normalization = sin(sqrtz*R0) / (sin(sqrtzShift*R0)+1i*tl*exp(1i*sqrtz*R0));
t_mat_ref(2, 1) = tl * normalization;

for j = 2 : numAtom
    t_mat_ref(:, j) = t_mat_ref(:, 1);
end
T_ref = diag(reshape(t_mat_ref, [ ], 1));

G_ref = (eye(2*numAtom) - G0*T_ref)^(-1) * G0;

% Check the decay rate of G_ref and G0
if (0)
x = zeros(numAtom, numAtom);
y_ref = zeros(numAtom, numAtom);
y0 = zeros(numAtom, numAtom);
for j = 1 : numAtom
    for n = 1 : numAtom
        x(j, n) = abs(j-n);
        y_ref(j, n) = norm(G_ref(2*j-1:2*j, 2*n-1:2*n), 'fro');
        y0(j, n) = norm(G0(2*j-1:2*j, 2*n-1:2*n), 'fro');
    end
end
x = reshape(x, [ ], 1);
y_ref = reshape(y_ref, [ ], 1);
y0 = reshape(y0, [ ], 1);
figure
semilogy(x, y_ref, 'o', 'color', [0.3010 0.7450 0.9330], ...
                    'linewidth', 2, 'markersize', 8, 'markerfacecolor', [0.3010 0.7450 0.9330])
title('$G_{\rm{ref}}$ off-diagonal decay', 'interpreter', 'latex', 'fontsize', 20)
figure
semilogy(x, y0, 'o', 'color', [0.3010 0.7450 0.9330], ...
                    'linewidth', 2, 'markersize', 8, 'markerfacecolor', [0.3010 0.7450 0.9330])
title('$G_{0}$ off-diagonal decay', 'interpreter', 'latex', 'fontsize', 20)
end


if (0)
x = zeros(numAtom, 1);
y_ref = zeros(numAtom, 1);
y0 = zeros(numAtom, 1);
for j = 1 : numAtom
    x(j) = abs(j-1);
    y_ref(j) = norm(G_ref(2*j-1:2*j, 1:2), 'fro');
    y0(j) = norm(G0(2*j-1:2*j, 1:2), 'fro');
end
save decay_data3.mat x y_ref y0
end


end