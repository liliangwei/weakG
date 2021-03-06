% 1D Laplace equation mixed finite element method (weak Galerkin)
% with Lagrangian multeplier
% Starson
% Sep 8

close all
clear all
format short e

x0 = 0; x1 = 1;   % start and end point
L = x1 - x0;      % total length
N = 10;          % number of interval
dx = L/N;         % interval size
dof = 2*N+1;      % degree of freedoms
ap = 10.;         % applied force on the RHS end

% test case 1 condition
% u0 = 0; uend=0;
% Decleration
u = zeros(1,dof); A = zeros(dof,dof); F = zeros(1,dof);
isBdary = false(dof,1); x = linspace(0,1,N+1);
T = [1:N]';

% Assemble global matrix
% Element stiffness matrix and load vector
Ak = dx/6*[2,1;1,2];
% Assemble global matrix
for i=1:N
    for m=1:2
        for n=1:2
   A(m+i-1,n+i-1) = A(m+i-1,n+i-1) + Ak(m,n);
        end
    end
    A(i,N+i+1)=1;A(i+1,N+i+1)=-1;
    A(N+i+1,i)=1;A(N+i+1,i+1)=-1;
end
A = -A;
% RHS
F(dof) = ap; F = F';
% Apply Boundary condition
isBdary(N+1) = true;
free = find(~isBdary);
% Calculate results
u(free) = A(free,free)\F(free);
% Error analysis
figure
plot(x,u(N+1:dof))
