% Follow FETI example then parrallel 9 elements poission equation
% with 3 individual subdomains. The center subdomain is unsupported
% Dirichlet boundary on each ending point

close all
clear all
format short e

L = 1;
N = 9;
dx = L/N;
x = linspace(0,1,N+1);
f1 = zeros(1,3); f3 = zeros(1,3); f2 = zeros(1,4);
for i=1:3
    
%Fk = [1/dx*(sin(2*pi*b)-sin(2*pi*a))-2*pi*cos(2*pi*a);
 %     1/dx*(sin(2*pi*a)-sin(2*pi*b))+2*pi*cos(2*pi*b)];
end
% unknowns of the subdomains
K1 = 1/dx*[2,-1,0;-1,2,-1;0,-1,1];
K2 = 1/dx*[1,-1,0,0;-1,2,-1,0;0,-1,2,0;0,0,-1,1];
K3 = 1/dx*[1,-1,0;-1,2,-1;0,-1,2];

% LU Gauss elimination  
[K1L,K1U] = lu(K1);
[K2L,K2U] = lu(K2);
[K3L,K3U] = lu(K3);

% Rigid body motion matrices
R1 = zeros(3,1); R3 = zeros(3,1); R2 = ones(4,1);

% Pseudo-inverse matrix of the stiffness matrices of particular domain.
K1P = inv(K1); K3P = inv(K3);
K2P = zeros(4,4);
K2P(1:3,1:3) = inv(K2(1:3,1:3));

% Connection matrices
epsilon1 = [0,0,-1;0,0,0];
epsilon2 = [1,0,0,0;0,0,0,-1];
epsilon3 = [0,0,0;1,0,0];

% Matrix G 
G = -epsilon1*R1-epsilon2*R2-epsilon3*R3;

% Vector e 