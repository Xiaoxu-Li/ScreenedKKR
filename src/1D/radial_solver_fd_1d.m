 function u = radial_solver_fd_1d(l, z, R, V0, N)
% Using finite difference method to solve the 1d schrodinger equation
% (-d^2/dr^2+V0(r)-z) u_l(r)=0 with boundary conditions 
% u_l(0)=0 and u_l(R)=1 for odd component (l=-1) and
% u'_l(0)=0 and u_l(R)=1 for even component (l=1)
%
% Note that this script is for symmetric potential V0
%
% 31/MAR/2023

h = R/N;
xdata = (0:N)'*h;
diag_ele = 2/h^2 - z + V0(xdata);
% construct the stiff matrix
H = diag(diag_ele) - diag(ones(N,1)/h^2, -1) - diag(ones(N,1)/h^2, 1);
H(1, :) = zeros(1, N+1);
if l == -1
    H(1,1) = 1;
else
    H(1,1) = 1;
    H(1,2) = -1;
end
H(end, :) = zeros(1, N+1);
H(end,end) = 1;

f = zeros(N+1, 1);
f(end) = 1;
% fprintf('rcond(H) is %e \n', rcond(H))

u = H\f;

end