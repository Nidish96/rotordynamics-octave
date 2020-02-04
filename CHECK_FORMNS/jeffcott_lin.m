clc
clear all 
addpath('../ROUTINES')

				% System Setup
%% Parameters
E   = 210e9;
rho = 7800;
g   = 9.81;

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

				% Linear Forced Response: Stepped SHBM - Linsolves
Nw = 2000;
Ws = linspace(Wcrit(1)/1000, Wcrit(5)*2, Nw);

Nd = size(Mb,1);
Ub = zeros(4, 3, Nw);
Upv = zeros(4, Nw);
X = zeros(Nd, 3, Nw);
for iw=1:Nw
  w = Ws(iw)*2*pi/60;
  E = HARMONICSTIFFNESS(Mb, Cb-w*Gb, Kb, w, [0; 1]);
  F = [Fgb; FCb*w^2; FSb*w^2];

  X(:, :, iw) = reshape(E\F, Nd, 3);

  Ub(:, :, iw) = Rb*X(:, :, iw);
  Upv(:, iw) = sqrt(Ub(:,1,iw).^2 + 0.5*sum(Ub(:,2:3,iw).^2,2));
end

		   % Linear Forced Response: Continued MHBM - Posed as NLsolves
Wstart = Wcrit(1)/1000*2*pi/60;
Wend = Wcrit(5)*2*2*pi/60;

ds = (Wend-Wstart)/2000;
Nt = 128;
h = [0 1];
Nhc = sum(h==0) + 2*sum(h~=0);

Copt = struct('Nmax', 1000, 'Display', 1);
Ubws = CONTINUE(@(Uw) ROTOR_RESFUN(Uw, Mb, Cb, Gb, Kb, Fgb, FCb, FSb, h, Nt), zeros(Nd*Nhc, 1), ...
		Wstart, Wend, ds, Copt);

C.Ws = Ubws(Nd*Nhc+1, :);
C.Uh = reshape(Ubws(1:end-1, :), Nd, Nhc, []);
C.Udyn_pv = Rb*sqrt(squeeze(sum(C.Uh(:, 2:end, :).^2, 2)));

%% Plotting
lw = 3;
%% Complare Continued Response
figure(100)
clf()
semilogy(C.Ws*60/2/pi, sqrt(0.5)*C.Udyn_pv(1, :), 'go-', 'LineWidth', 1); hold on
semilogy(Ws, Upv(1, :), 'b.', 'LineWidth', lw); hold on

%% Response
figure(1)  % Mid point RMS
clf()

subplot(2, 1, 1)
plot(Ws, Upv(1, :), 'b-', 'LineWidth', lw); hold on
plot(Ws, Upv(2, :), 'r-', 'LineWidth', lw); hold on
set(gca, 'YScale', 'log');
legend('Horizontal Deflection', 'Vertical Deflection')

for ic=[1 3 5]
  plot([Wcrit(ic) Wcrit(ic)], [1e-24 1e-3], 'k--', 'LineWidth', lw)
end
xlim([Ws(1) Ws(end)])
ylim([1e-8 1e-3])
xlabel('Spin Speed (rpm)');
ylabel('RMS Amplitude (m)');

subplot(2, 2, 3)
semilogy(Ws, abs(Ub(1, 1, :)), 'k-', 'LineWidth', lw); hold on
semilogy(Ws, sqrt(sum(Ub(1, 2:3, :).^2,2)), 'b-', 'LineWidth', lw)
title('Horizontal Deflection')

legend('Static', 'Harmonic')
xlim([Ws(1) Ws(end)])
ylim([1e-8 1e-3])
xlabel('Spin Speed (rpm)');
ylabel('Deflection Amplitude (m)')

subplot(2, 2, 4)
semilogy(Ws, abs(Ub(2, 1, :)), 'k-', 'LineWidth', lw); hold on
semilogy(Ws, sqrt(sum(Ub(2, 2:3, :).^2,2)), 'r-', 'LineWidth', lw)
title('Vertical Deflection')

legend('Static', 'Harmonic')
xlim([Ws(1) Ws(end)])
ylim([1e-8 1e-3])
xlabel('Spin Speed (rpm)');
ylabel('Deflection Amplitude (m)')

%% Plot Deflection Shape
Nx = 2;
Ny = 3;
IWs = fix(linspace(1, Nw, Nx*Ny));

figure(2)
clf()

for i=1:Nx*Ny
  iw = IWs(i);
  % sc = 2/max(max(abs(X(:, 2:3, iw))));
  sc = 7e3;

  subplot(Nx, Ny, i)
  DEPICTROTOR(Xs, Ris, Ros, Disc, Lb*X(:, :, iw)*sc);
  axis equal
  title(sprintf('Spin Speed: %.2f rpm', Ws(iw)))
end
