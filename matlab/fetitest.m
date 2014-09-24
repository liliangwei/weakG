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
f1 = zeros(1,4); f3 = zeros(1,4); f2 = zeros(1,4);
for i=1:N
    a = x(i); b = x(i+1);   
    Fk = [1/dx*(sin(2*pi*b)-sin(2*pi*a))-2*pi*cos(2*pi*a);
      1/dx*(sin(2*pi*a)-sin(2*pi*b))+2*pi*cos(2*pi*b)];
  if(i<4)
      for m=1:2
          f1(m+i-1) = f1(m+i-1) + Fk(m);
      end
  elseif(i<7)
      for m=1:2
          f2(m+i-4) = f2(m+i-4) + Fk(m);
      end
  else
      for m=1:2
          f3(m+i-7) = f3(m+i-7) + Fk(m);
      end
  end
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
G = [-epsilon2*R2];

% Vector e 
e = -dot(R2,f2);

% lambda0
lambda0 = G*inv(dot(G',G))*e;

% F matrix
F = epsilon1*K1P*epsilon1' + epsilon2*K2P*epsilon2' + epsilon3*K3P*epsilon3';

% Vector g
g = epsilon1*K1P*f1(2:4)' + epsilon2*K2P*f2' + epsilon3*K3P*f3(1:3)';

% simple one-dimensional example leads to the specially reduced system of
% equtions where no iteration is necessary
alpha = inv(dot(G',G))*G'*(g-F*lambda0);

% Displacemenets of the subdomains are obtained from two contributions,
% dker and dinfi.
% The dinfi are obtained with the help of pseudo-inverse matrices

dinfi1 = K1P*(f1(2:4)' - epsilon1' * lambda0);
dinfi2 = K2P*(f2' - epsilon2' * lambda0);
dinfi3 = K3P*(f3(1:3)' - epsilon3' * lambda0);

% Contribuions from rigid body motions 
dker1 = R1*0;
dker2 = R2*alpha;
dker3 = R3*0;

% The final subdomain displacement vectors are
d1 = dinfi1 + dker1;
d2 = dinfi2 + dker2;
d3 = dinfi3 + dker3;

% plot 
u1 = zeros(1,4); u2 = zeros(1,4); u3 = zeros(1,4);
u1(2:4) = d1; u2 = d2; u3(1:3) = d3;

figure
plot(x(1:4),u1);
hold on
plot(x(4:7),u2,'r-');
hold on
plot(x(7:10),u3,'g-');
