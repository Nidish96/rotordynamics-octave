clc
clear all
addpath('../ROUTINES')

%% Parameters
E = 2.078e11;
rho = 7806;

Ne = 9;
Nn = Ne+1;

L = 1.0;
D = 1e-2;

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
    [Me, Ge, Ke, ~, MTe, MRe] = ROTELMATS(E, rho, L/Ne, 0, D/2, 0);
    
    is = (e-1)*4+(1:8);
    M(is, is) = M(is, is) + Me;
    G(is, is) = G(is, is) + Ge;    
    K(is, is) = K(is, is) + Ke;
end

%% Boundary Conditions
Lb = eye(size(M));
Lb(:, 1:4) = [];

Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;
Gb = Lb'*G*Lb;

%% Verify
A = pi*D^2/4;
I = pi*D^4/64;
mmlp = sqrt(E*I/(rho*A*L^4));

Ws = sort(sqrt(eig(Kb, Mb)));

sqrt(Ws(1:10)/mmlp)
 
