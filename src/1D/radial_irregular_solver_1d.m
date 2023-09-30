 function u = radial_irregular_solver_1d(l, z, R0, V0, N)
% Using finite difference method to solve an irregular solution of 
% the radial schrodinger equation 
% (-d^2/dr^2+V(r)-z)u_l(r)=0 with boundary conditions 
% u_l(R0)=h_l(sqrt{z}R0) and u'_l(R0)=sqrt{z} h'_l(sqrt{z}R0)
%
% 31/MAR/2023

h = R0/N;
xdata = (0:N)'*h;
diag_ele = 2/h^2 - z + V0(xdata);
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

if l == 1
    f(N) = 1i*sqrtz*exp(1i*sqrtz*R0)*h;
    f(N+1) = exp(1i*sqrtz*R0);
else
    f(N) = sqrtz*exp(1i*sqrtz*R0)*h;
    f(N+1) = -1i*exp(1i*sqrtz*R0);
end

u = H\f;
% fprintf('rcond(H) is %e \n', rcond(H))


end