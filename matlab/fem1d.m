% 1D Poisson equation finite element method
% Standard Galerkin method
% Starson
% Sep 6

close all
clear all
format short e

x0 = 0; x1 = 1;   % start and end point
L = x1 - x0;      % total length
N = 5;            % number of interval
dx = L/N;         % interval size

% test case 1 condition
% u0 = 0; du1/dt = 0; f=1;
% Decleration
u = zeros(1,N+1); A = zeros(N+1,N+1); F = zeros(1,N+1);
isBdary = false(N+1,1); x = linspace(0,1,N+1);
% Element stiffness matrix and load vector
Ak = 1/dx*[1,-1;-1,1];
Fk = dx/2*[1;1];
% Assemble global matrix
for i=1:N
    for m=1:2
        for n=1:2
   A(m+i-1,n+i-1) = A(m+i-1,n+i-1) + Ak(m,n);
        end
        F(m+i-1) = F(m+i-1) + Fk(m);
    end
end
F = -2*F';

% Apply Boundary condition
isBdary(1) = true;isBdary(N+1) = true;
free = find(~isBdary);

% Calculate results
u(free) = A(free,free)\F(free);
u_exact = -x.*(1-x);

% Error analysis
error = sqrt(sum((u_exact-u).^2))/N;
disp(error);

figure
plot(x,u);
hold on
plot(x,u_exact,'r-');
hold off