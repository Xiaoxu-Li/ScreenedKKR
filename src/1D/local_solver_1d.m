function [tl, A, B] = local_solver_1d(l, z, R, V0, N)
% Solves (-d^2/dx^2 + V(x)) u(x) = Eu(x) locally in a single cell
% Construct the t-matrix and two integral values A and B
% A = \int_{0}^{R} |u|^2 and B = \int_{0}^{R} u'v
% where Omega is the current cell, and u and v are 
% regular and irregular local radial solutions respectively 
% with the angular momentum l for a muffin-tin system
% 
% We first solve a radial schrondinger equation using FDM
% and then obtain the t-matrix via the integral definition
% At last, we need to normalize the t-matrix, local solutions to match the BC
%
% 05/APR/2023

u = radial_solver_fd_1d(l, z, R, V0, N);
v = radial_irregular_solver_1d(l, z, R, V0, N);

h = R/N; 
xgrid = (1:N)'*h;
u = u(2:end);
v = v(2:end);

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
    normalization = cos(sqrtz*R) / (1-tl*exp(1i*sqrtz*R));
else
    normalization = sin(sqrtz*R) / (1+1i*tl*exp(1i*sqrtz*R));
end
tl = tl * normalization;

u = u * normalization;
A = u'*u*h;
B = u'*v*h;

end