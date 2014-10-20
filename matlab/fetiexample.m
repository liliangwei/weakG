% Follow FETI example then parrallel linear elastic problem
% with 3 individual subdomains.
% Dirichlet boundary on LHS ending.

close all
clear all
format short e

l = 1.0;                  % length
N = 12;                   % number of elements
x = linspace(0,l,N+1);    % x coordinate
A = 0.1;                  % cross-sectional area
E = 0.1;                  % Young's modulu of rubber
F = 1.0;                  % applied force on RHS end

% First part is regular single domain linear elastic problem

k = E*A/l*[1,-1;-1,1];    % element stiffness matrix
f_s = zeros(N+1,1); f_s(N+1) = F;  % applied force 
isBdary = false(N+1,1); isBdary(1) = true;  % boundary condition
k_s = zeros(N+1,N+1);        % global stiffness matrix
d_s = zeros(N+1,1);          % displacement
for i=1:N
    for m=1:2
        for n=1:2
            k_s(m+i-1,n+i-1) = k_s(m+i-1,n+i-1) + k(m,n);
        end
    end
end

free = find(~isBdary);
d_s(free) = k_s(free,free)\f_s(free); % calculate the displacement

figure
plot(x,d_s,'-bs','LineWidth',2);
title('Displacement distribution for single domain');

% Second part is multi-subdomain 

n_domain = 3;   % number of subdomains
N_sub = (N+n_domain)/n_domain;  % dof for each subdomain
d1 = zeros(N_sub,1); d2 = zeros(N_sub,1); d3 = zeros(N_sub,1);
k_sub = zeros(N_sub,N_sub);
f1 = zeros(N_sub-1,1); f2 = zeros(N_sub,1); f3=f2; f3(N_sub) = F;
for i=1:N_sub-1
    for m=1:2
        for n=1:2
            k_sub(m+i-1,n+i-1) = k_sub(m+i-1,n+i-1) + k(m,n);
        end
    end
end
% stiffness matrix for each subdomain
k1 = k_sub(2:N_sub,2:N_sub); k2 = k_sub; k3 = k_sub;

% rigid body motion vectors
R1 = zeros(N_sub-1,1); R2 = ones(N_sub,1); R3 = ones(N_sub,1);

% pseudo-inverse matrices
K2P = zeros(N_sub,N_sub);
K1P = inv(k1); K2P(1:N_sub-1,1:N_sub-1) = inv(k2(1:N_sub-1,1:N_sub-1));
K3P = K2P;

% assembling matrices
epsilon1 = [0,0,0,-1;0,0,0,0];
epsilon2 = [1,0,0,0,0;0,0,0,0,-1];
epsilon3 = [0,0,0,0,0;1,0,0,0,0];

% matrix G
G = [-epsilon2*R2,-epsilon3*R3];

% vector e
e = [(-R2'*f2),(-R3'*f3)];
e = e';

% initial approximation of the vector of Lagrange multipliers
lambda0 = G*inv(G'*G)*e;

% The F matrix
F_sub = epsilon1*K1P*epsilon1'+epsilon2*K2P*epsilon2'+epsilon3*K3P*epsilon3';

g = zeros(2,1);

% this special simple one-dimensional example leads to the specially
% reduced system of equations where no iteration is necessary.

alpha = inv(G'*G)*G'*(g-F_sub*lambda0);

% displacement obtained with the pseudo-inverse matrices
dinfi1 = K1P*(f1-epsilon1'*lambda0);
dinfi2 = K2P*(f2-epsilon2'*lambda0);
dinfi3 = K3P*(f3-epsilon3'*lambda0);

% displacement contributions from rigid body motions
dker1 = R1;
dker2 = R2*alpha(1);
dker3 = R3*alpha(2);

d1(2:N_sub) = dinfi1+dker1;
d2 = dinfi2+dker2;
d3 = dinfi3+dker3;

figure
plot(x(1:5),d1,'-rs','LineWidth',3);
hold on
plot(x(5:9),d2,'--gs','LineWidth',2);
hold on
plot(x(9:13),d3,'-.ms','LineWidth',1);
hold off
title('Displacement distribution for multiple subdomains')

