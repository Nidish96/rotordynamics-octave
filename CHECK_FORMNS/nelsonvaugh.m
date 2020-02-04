clc
clear all
addpath('../ROUTINES')

%% Parameters
E = 2.078e11;
rho = 7806;

Xs = [0.00, 1.27, 5.08, 7.62, 8.89, 10.16, 10.67, 11.43, 12.70, 13.46, ...
      16.51, 19.05, 22.86, 26.67, 28.70, 30.48, 31.50, 34.54, 35.81]*1e-2;
Ris = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.52, 1.78, 0.00, 0.00, 0.00, ...
       0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.52]*1e-2;
Ros = [0.51, 1.02, 0.76, 2.03, 2.03, 3.30, 3.30, 2.54, 2.54, 1.27, 1.27, ...
       1.52, 1.52, 1.27, 1.27, 3.81, 2.03, 2.03]*1e-2;

Nn = length(Xs);
Ne = Nn-1;

g = 9.81; % gravity

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
  [Me, Ge, Ke, ~] = ROTELMATS(E, rho, Xs(e+1)-Xs(e), Ris(e), Ros(e), g);

  is = (e-1)*4+(1:8);
  M(is, is) = M(is, is) + Me;
  G(is, is) = G(is, is) + Ge;
  K(is, is) = K(is, is) + Ke;
end

%% Discrete Discs
Disc = struct('m', 1.401, 'J', 0.00136, 'Jz', 0.00203, 'ax', 0.635e-3);
M((5-1)*4+(1:4), (5-1)*4+(1:4)) = M((5-1)*4+(1:4), (5-1)*4+(1:4)) + ...
				  Disc.m*diag([1 1 0 0]) + ...
				  Disc.J*diag([0 0 1 1]);

G((5-1)*4+(3:4), (5-1)*4+(3:4)) = G((5-1)*4+(3:4), (5-1)*4+(3:4)) + ...
				  Disc.Jz*[0 -1;1 0];

FC((5-1)*4+(1:4)) = FC((5-1)*4+(1:4)) + [Disc.m*Disc.ax; 0; 0; 0];
FS((5-1)*4+(1:4)) = FS((5-1)*4+(1:4)) + [0; Disc.m*Disc.ax; 0; 0];
				  
%% Bearings
Ka = zeros(size(K));
Kb = zeros(size(K));

Ka((11-1)*4+(1:2), (11-1)*4+(1:2)) = 4.378e7*eye(2);
Ka((15-1)*4+(1:2), (15-1)*4+(1:2)) = 4.378e7*eye(2);

Kb((11-1)*4+(1:2), (11-1)*4+(1:2)) = [3.503e7 -0.8756e7;-0.8756e7 3.503e7];
Kb((15-1)*4+(1:2), (15-1)*4+(1:2)) = [3.503e7 -0.8756e7;-0.8756e7 3.503e7];

Ka = K + Ka;
Kb = K + Kb;

%% QEP
W = 1000*2*pi/60;

[V, Z] = polyeig(Kb, -G*W, M);
[~, si] = sort(abs(Z));
Z = Z(si);
V = V(:, si);

reshape(imag(Z(1:2:12))*60/2/pi, 2, [])'

figure(1)
clf()
DEPICTROTOR(Xs, Ris, Ros, []);
axis equal
