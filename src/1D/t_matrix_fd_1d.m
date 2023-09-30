function tl = t_matrix_fd_1d(l, z, R0, V0, N)
% Consider (-d^2/dx^2 + V(x)) u(x) = Eu(x) within 1D MST
% Construct the t-matrix for a muffin-tin system
% 
% We first solve a radial schrondinger equation using FDM,
% and then obtain the t-matrix via the integral definition
% At last, we need to normalize the t-matrix to match the BC
%
% 29/MAR/2023

u = radial_solver_fd_1d(l, z, R0, V0, N);

h = R0/N; 
xgrid = (1:N)'*h;
u = u(2:end);

% Guarantee the imaginary part of sqrt(z) is positive
if imag(sqrt(z))>=0
    sqrtz = sqrt(z); 
else
    sqrtz = -sqrt(z);
end

if l == 1
    bessel_1d = cos(sqrtz*xgrid);
else
    bessel_1d = sin(sqrtz*xgrid);
end
tl = sum( bessel_1d .* V0(xgrid) .* u) * h;
tl = tl / (1i * sqrtz);

if l == 1
    normalization = cos(sqrtz*R0) / (1-tl*exp(1i*sqrtz*R0));
else
    normalization = sin(sqrtz*R0) / (1+1i*tl*exp(1i*sqrtz*R0));
end
tl = tl * normalization;

end