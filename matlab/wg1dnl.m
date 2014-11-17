% 1D linear elastic equation weak Galerkin finite element method
% Starson


close all
clear all
format short e

x0 = 0; x1 = 1;   % start and end point
L = x1 - x0;      % total length
N = 5;          % number of interval
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
Moo = Zk'*Ak*Zk + Zk'*Ak*Zk;
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

% ------------------------------------------------------------
% A new implementation of element matrix

newM = zeros(3,3);
newM(1,1) = Mbb(1,1); newM(1,3) = Mbb(1,2); 
newM(3,1) = Mbb(2,1); newM(3,3) = Mbb(2,2);
newM(2,1) = Mbo(1,1); newM(3,2) = Mbo(2,1);
newM(1,2) = Mob(1,1); newM(2,3) = Mbo(2,1);
newM(2,2) = Moo(1,1);

newA = zeros(dof,dof);
for i=1:N
    for p=1:3
        for q=1:3
            newA(p+(i-1)*2,q+(i-1)*2) = newA(p+(i-1)*2,q+(i-1)*2) + newM(p,q);
        end
    end
end

newu = zeros(1,dof); isBdary2 = false(dof,1);
isBdary2(1) = true;
nfree = find(~isBdary2);
newu(nfree) = newA(nfree,nfree)\F(nfree);

figure
plot(x,u(N+1:dof));
grid on