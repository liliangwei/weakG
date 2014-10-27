% 1D linear elastic equation weak Galerkin finite element method
% Starson


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
Dk = dx/6*[2,-1;-1,2]; 
Zk = [1;1]; Tk = [1,0;0,1];
Ak = 2/dx*[2,1;1,2];
% Assemble global matrix
Moo = Zk'*Ak*Zk;
Mob = -Zk'*Ak*Tk; Mbo = Mob';
Mbb = Tk'*Ak*Tk;
A = A+sparse(T,T,Moo,dof,dof);
Iu = [T,T]; Iv = [N+T,N+T+1];
val = zeros(N,2);
val(:,1) = Mob(1,1);val(:,2) = Mob(1,2);

A = A+sparse(Iu,Iv,val,dof,dof);
A = A+sparse(Iv,Iu,val',dof,dof);
for i=N+1:dof-1
    for m=1:2
        for n=1:2
   A(m+i-1,n+i-1) = A(m+i-1,n+i-1) + Mbb(m,n);
        end
    end
end

% RHS
F(dof) = ap; F = F';
% Apply Boundary condition
isBdary(N+1) = true;
free = find(~isBdary);

% Calculate results
u(free) = A(free,free)\F(free);

figure
plot(x,u(N+1:dof));
grid on