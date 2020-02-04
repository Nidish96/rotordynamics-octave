clc
clear all
addpath('../ROUTINES')

%% Parameters
E = 2.078e11;
rho = 7806;

Ne = 5;
Nn = Ne+1;

L = 127e-2;  % 127 cm shaft

Ri = 0;
Ro = 10.16e-2/2;  % 10.16 cm diameter

g = 0;  % gravity

%% Matrices
M = zeros(Nn*4);
K = zeros(Nn*4);
G = zeros(Nn*4);
FC = zeros(Nn*4, 1);
FS = zeros(Nn*4, 1);

Me = zeros(8);
Ke = zeros(8);
Ge = zeros(8);
for e=1:Ne
  [Me, Ge, Ke, ~] = ROTELMATS(E, rho, L/Ne, Ri, Ro, g);
    
  is = (e-1)*4+(1:8);
  M(is, is) = M(is, is) + Me;
  G(is, is) = G(is, is) + Ge;    
  K(is, is) = K(is, is) + Ke;
end

%% Boundary Nodes: Linear Springs
bn_ids = [1 Nn];
bn_ddofs = reshape((bn_ids-1)*4+[1;2],[],1);

K(bn_ddofs, bn_ddofs) = K(bn_ddofs, bn_ddofs) + eye(4)*1.753e7;

%% QEP
W = 4000*2*pi/60;

[V, Z] = polyeig(K, -G*W, M);
[~, si] = sort(abs(Z));
Z = Z(si);
V = V(:, si);

imag(Z(1:2:16))
