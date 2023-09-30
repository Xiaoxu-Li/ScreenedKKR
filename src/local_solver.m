function [tl, A, B] = local_solver(l, z, R0, V0, N)
% Solves (-\nabla + V(r)) u(r) = Eu(r) locally in a single cell
% Construct the t-matrix and two integral values A and B
% A = \int_{0}^{R0} r^2 |u|^2 and B = -i\sqrt(z) \int_{0}^{R0} r^2 u'v
% where Omega is the current cell and u and v are 
% regular and irregular local radial solutions respectively 
% with the angular momentum l for a muffin-tin system
% 
% We first solve a radial schrondinger equation using FDM
% and then obtain the t-matrix via the integral definition
% At last, we need to normalize the t-matrix, local solutions to match the BC
%
% 28/JAN/2022

u = radial_solver_fd(l, z, R0, V0, N);
v = radial_irregular_solver(l, z, R0, V0, N);
% Note that u and v corresponding to rR_l(r)

h = R0/N; 
xgrid = (1:N)'*h;
u = u(2:end);
% Note that j0(x)=sinx/x needs to be redefined at the origin, 
% thus we use the right endpoints
v = v(2:end);

% Guarantee the imaginary part of sqrt(z) is positive
if imag(sqrt(z))>=0
    sqrtz = sqrt(z); 
else
    sqrtz = -sqrt(z);
end
tl = sum( spherical_bessel(l, sqrtz*xgrid) .* V0(xgrid) .* xgrid .* u) * h;
tl = -1i * sqrtz * tl;

normalization = spherical_bessel(l, sqrtz*R0) ...
    / (1-tl*spherical_hankel(l, sqrtz*R0));
tl = tl * normalization;

% Normalization for the local regular solutions
u_normalization = spherical_bessel(l, sqrtz*R0) + tl * spherical_hankel(l, sqrtz*R0);
u = u * u_normalization;
A = u'*u*h;
B = 1i*sqrtz*u'*v*h;

end