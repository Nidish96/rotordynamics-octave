clc
clear all 
addpath('../ROUTINES')

analyze = true;  % false for load data

				% System Setup
%% Parameters
E   = 210e9;
rho = 7800;
g   = 9.81*0;

Xs = linspace(0, 127e-3, 5);
Ris = zeros(1, 4);
Ros = 10.16e-3/2*ones(1, 4);

Disc = DISC(3, Ros(3), 25e-3, 5e-3, rho, 6e-6, 0, g);

Nn = length(Xs);
Ne = Nn-1;

%% Matrices
M = zeros(Nn*4);
G = zeros(Nn*4);
K = zeros(Nn*4);
Fg = zeros(Nn*4, 1);
FC = zeros(Nn*4, 1);
FS = zeros(Nn*4, 1);

Me = zeros(8);
Ge = zeros(8);
Ke = zeros(8);
Fge = zeros(8, 1);
for e=1:Ne
  [Me, Ge, Ke, Fge] = ROTELMATS(E, rho, Xs(e+1)-Xs(e), Ris(e), Ros(e), g);

  is = (e-1)*4+(1:8);
  M(is, is) = M(is, is) + Me;
  G(is, is) = G(is, is) + Ge;
  K(is, is) = K(is, is) + Ke;
  Fg(is) = Fg(is) + Fge;
end

%% Disc
[Md, Gd, FCd, FSd, Fgd] = Disc.MATS();
dis = (Disc.nd-1)*4+(1:4);
M(dis, dis) = M(dis, dis) + Md;
G(dis, dis) = G(dis, dis) + Gd;
FC(dis) = FC(dis) + FCd;
FS(dis) = FS(dis) + FSd;
Fg(dis) = Fg(dis) + Fgd;

%% Boundary Conditions 
Lb = eye(Nn*4); Lb(:, [1:4 end-3:end]) = [];
Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;
Gb = Lb'*G*Lb;
FCb = Lb'*FC;  FSb = Lb'*FS;
Fgb = Lb'*Fg;

%% Extraction
R = eye(Nn*4); R = R((Disc.nd-1)*4+(1:4), :);
Rb = R*Lb;

%% Critical Speeds
[V, Ws] = eig(Kb, Mb);
[Ws, si] = sort(sqrt(diag(Ws)));  V = V(:, si);
Wcrit = sort(Ws*60/2/pi);

%% Linear Damping (Rayleigh Quotients)
Zt_req = [0.02e-2; 0.02e-2; 0.04e-2; 0.04e-2];
PHI = [0.5./Ws(1:length(Zt_req)) Ws(1:length(Zt_req))];  Y = Zt_req;
ab = PHI\Y;
Cb = ab(1)*Mb + ab(2)*Kb;

%% Clearance Nonlinearities
nl_els = struct('type', 'clearance', 'pars', [1e-6; 1e8], ...
		'sel_shape', Lb((Disc.nd-1)*4+(1:2), :), ...
		'f_shape', Lb((Disc.nd-1)*4+(1:2), :)');

%% Harmonic Balance Continuation
h   = [0 1 3 5 7];
Nd  = size(Mb, 1);
Nhc = sum(h==0) + 2*sum(h~=0);
Nt  = 128;

Wstart = Wcrit(1)/1000*2*pi/60;
Wend = Wcrit(5)*2*2*pi/60;

ds = (Wend-Wstart)/2000;

Copt = struct('Nmax', 1000, 'Display', 1, 'angopt', 1e-6, 'dsmin', 10);

if analyze
  E0 = HARMONICSTIFFNESS(Mb, Cb-Wstart*Gb, Kb, Wstart, h);
  U0 = E0\([Fgb; Wstart^2*[FCb; FSb]; zeros(Nd*(Nhc-3),1)]);

  Ubws = CONTINUE(@(Uw) ROTOR_RESFUN(Uw, Mb, Cb, Gb, Kb, Fgb, FCb, FSb, h, Nt, nl_els), ...
		  zeros(Nd*Nhc, 1), Wstart, Wend, ds, Copt);
  
  C.Ws = Ubws(Nd*Nhc+1, :);
  C.Uh = reshape(Ubws(1:end-1, :), Nd, Nhc, []);
  C.Udyn_pv = Rb*sqrt(squeeze(C.Uh(:, 1, :).^2+0.5*sum(C.Uh(:, 2:end, :).^2, 2)));

  save('./DATS/jeffcott_clearance_nlin.mat', 'Ubws', 'C', 'nl_els', 'h')
else
  load('./DATS/jeffcott_clearance_nlin.mat')
end

%% Plot
lw = 3;
figure(1)
clf()
semilogy(C.Ws, C.Udyn_pv(1, :), 'b-', 'Linewidth', lw); hold on
xlim([0 1.5e5])
# semilogy(C.Ws, Rb(1,:)*abs(squeeze(0.5*sum(C.Uh(:, 6:7, :).^2,2))), 'k--')
