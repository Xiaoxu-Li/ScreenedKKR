 function u = radial_irregular_solver(l, z, R0, V0, N)
% Using finite difference method to solve an irregular solution of 
% the radial schrodinger equation
% (-d^2/dr^2+l(l+1)/r^2+V(r)-z)u_l(r)=0 where u_l(r)=rR_l(r)
% with boundary conditions u_l(R0)=R0*h_l(sqrt{z}R0) 
% and u'_l(R0)=R0*sqrt{z} h'_l(sqrt{z}R0)
%
% Note that the radial equation corresponds to (-\nabla + V(r)) u(r) = zu(r)
%
% 03/FEB/2023

h = R0/N;
xdata = (0:N)'*h;
diag_ele = 2/h^2 - z + l*(l+1)./(xdata.^2) + V0(xdata);
% construct the stiff matrix
H_temp = diag(diag_ele) - diag(ones(N,1)/h^2, -1) - diag(ones(N,1)/h^2, 1);
H = zeros(N+1, N+1);
H(1:N-1, :) = H_temp(2:N, :);
H(end-1, end-1) = -1;
H(end-1,end) = 1;
H(end,end) = 1;

% Guarantee the imaginary part of sqrt(z) is positive
if imag(sqrt(z))>=0
    sqrtz = sqrt(z); 
else
    sqrtz = -sqrt(z);
end
f = zeros(N+1, 1);
f(N) = R0 * (spherical_hankel(l, sqrtz*R0) - spherical_hankel(l, sqrtz*(R0-h)));
f(N+1) = R0 * spherical_hankel(l, sqrtz*R0);

u = H\f;
% fprintf('rcond(H) is %e \n', rcond(H))


end