% 1D Laplace equation weak Galerkin finite element method
% with Lagrangian multeplier
% Starson
% Sep 8

close all
clear all
format short e

x0 = 0; x1 = 1;   % start and end point
L = x1 - x0;      % total length
N = 3;          % number of interval
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
% M00
A = A + sparse(T,T,12*dx,dof,dof);
% M0b & Mb0
Iu = [T,T];
Iv = N+[T,T+1];
val = dx*dx*[-6*ones(N,1),6*ones(N,1)];
A = A + sparse(Iu,Iv,val,dof,dof);
A = A + sparse(Iv,Iu,val',dof,dof);
% Mbb
Iu = N+[T,T;T+1,T+1];
Iv = Iu';
val = dx*dx*dx*[3*ones(N,1),-3*ones(N,1);-3*ones(N,1),3*ones(N,1)];
A = A + sparse(Iu,Iv,val,dof,dof);

% RHS
F(dof) = ap; F = F';
% Apply Boundary condition
isBdary(dof) = true;
free = find(~isBdary);
A(N+1,:) = 0.;A(N+1,N+1) = 1;
ai = inv(A);
% Calculate results
%u(free) = A(free,free)\F(free);

% Error analysis

