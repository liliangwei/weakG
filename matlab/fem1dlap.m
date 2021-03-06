% 1D Poisson equation finite element method
% Standard Galerkin method
% Starson


close all
clear all
format short e

x0 = 0; x1 = 1;   % start and end point
L = x1 - x0;      % total length
N = 5;            % number of interval
dx = L/N;         % interval size
ap = 10.;         % applied force on the RHS ending

% test case 1 condition
% u0 = 0; du1/dt = f; f=1;
% Decleration
u = zeros(1,N+1); A = zeros(N+1,N+1); F = zeros(1,N+1);
isBdary = false(N+1,1); x = linspace(0,1,N+1);
% Element stiffness matrix and load vector
Ak = 1/dx*[1,-1;-1,1];
% Assemble global matrix
for i=1:N
    for m=1:2
        for n=1:2
   A(m+i-1,n+i-1) = A(m+i-1,n+i-1) + Ak(m,n);
        end
        a = x(i); b = x(i+1);
    end
end
F(N+1) = ap;

% Apply Boundary condition
 isBdary(1) = true; isBdary(N+1) = true;
free = find(~isBdary);
A(1,:) = 0; A(1,1) = 1;
A(N+1,:) = 0; A(N+1,N+1) = 1;
F(1) = 5; F(free) = 8;

% Calculate results
%u(free) = A(free,free)\F(free)';
u = A \ F';
%u_exact = -sin(2*pi*x);

u(1) = 5;

% Error analysis
%error = sqrt(sum((u_exact-u).^2))/N;
%disp(error);

figure
plot(x,u);
hold on
%plot(x,u_exact,'r-');
hold off
grid on