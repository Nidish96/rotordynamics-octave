clc
clear all 
addpath('../ROUTINES')

analyze = true;
depict = false;
lissajous = true;
nlin = 'stator';
mu = 0.8;
disc = true;
upordown = 'up';
h = [0 1 2 3 4 5 6 7 8 9];
Nt = 2^10;
%% Parameters
E   = 190e9;
rho = 7.85e3;
g   = 0;

				% Model: [  =|===|=  ]
Lmid = 1697;
Les = [76.35 (Lmid/3+6.35) Lmid/3 (Lmid/3+3.175) 33.175]*1e-3;
Xs = [0 cumsum(Les)];
Nn = length(Xs);
Ne = length(Les);

Ris = zeros(size(Xs));
Ros = 0.25*25.4e-3+Ris;  % 1/2 inch diameter

Disc(1) = DISC(2, Ros(2), 3.5*25.4e-3, 0.5*25.4e-3, 7.8e3, 6e-3, 0, g);  % Cast iron flywheel
if disc
  Disc(2) = DISC(Nn-1, Ros(Nn-1), 2.5*25.4e-3, 0.25*25.4e-3, 8e3, 1e-3, 0, g);  % 303 SS large disc
  ## Disc(2) = DISC(Nn-1, Ros(Nn-1), 2*25.4e-3, 0.25*25.4e-3, 8e3, 0, 0, g);  % 303 SS small disc
end

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
for i=1:length(Disc)
  [Md, Gd, FCd, FSd, Fgd] = Disc(i).MATS();
  dis = (Disc(i).nd-1)*4+(1:4);
  M(dis, dis) = M(dis, dis) + Md;
  G(dis, dis) = G(dis, dis) + Gd;
  FC(dis) = FC(dis) + FCd;
  FS(dis) = FS(dis) + FSd;
  Fg(dis) = Fg(dis) + Fgd;
end

%% Boundary Conditions
Lb = eye(Nn*4);  Lb(:, 1:4) = [];  % Cantilevered
Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;
Gb = Lb'*G*Lb;
FCb = Lb'*FC;  FSb = Lb'*FS;
Fgb = Lb'*Fg;

%% Output Point
Rb = Lb(end-3:end,:);

%% Non-linearity
if nlin
  nl_els = struct('type', nlin, 'pars', [0.5*25.4e-3, 1e8, mu], ...
		  'sel_shape', Lb((Nn-1-1)*4+(1:2), :), ...
		  'f_shape', Lb((Nn-1-1)*4+(1:2), :)');
else
  nl_els = [];
end

%% Critical Speeds
[V, Ws] = eig(Kb, Mb);
[Ws, si] = sort(sqrt(diag(Ws)));  V = V(:, si);
Wcrit = sort(Ws*60/2/pi);
fprintf('%f\n', Wcrit(1))

%% Linear Damping (Rayleigh Quotients)
Zt_req = [0.5e-2; 0.5e-2; 0.4e-2; 0.4e-2];
PHI = [0.5./Ws(1:length(Zt_req)) Ws(1:length(Zt_req))];  Y = Zt_req;
ab = PHI\Y;
Cb = ab(1)*Mb + ab(2)*Kb;

	 % Linear Forced Response: Continued SHBM - Posed as NL Solves
if strcmp(upordown,'up')
  Wstart = 5*2*pi/60;  % 5 rpm
  Wend = 250*2*pi/60;  % 250 rpm
else
  Wstart = 250*2*pi/60;  % 5 rpm
  Wend = 5*2*pi/60;  % 250 rpm
end

Nd = size(Mb, 1);

ds = (Wend-Wstart)/100;
Nhc = sum(h==0) + 2*sum(h~=0);
Nhmax = max(h);

Copt = struct('Nmax', 1000, 'Display', 1, 'angopt', 1e-4);
if analyze
  U0 = HARMONICSTIFFNESS(Mb, Cb, Kb, Wstart, h)\[Fgb; FCb*Wstart^2; FSb*Wstart^2; zeros(Nd*(Nhc-3),1)];
  Ubws = CONTINUE(@(Uw) ROTOR_RESFUN(Uw, Mb, Cb, Gb, Kb, Fgb, FCb, FSb, h, Nt, nl_els), ...
		  U0, Wstart, Wend, ds, Copt);

  D1 = HARMONICSTIFFNESS(0, 1, 0, 1, h);
  C.Ws = Ubws(end, :);
  C.Uh = reshape(Ubws(1:end-1, :), Nd, Nhc, []);
  C.Vh = reshape(kron(D1, eye(Nd))*Ubws(1:end-1,:).*C.Ws, Nd, Nhc, []);
  C.Udyn_amp = Rb*sqrt(squeeze(sum(C.Uh(:, 2:end, :).^2, 2)));
  C.Vdyn_amp = Rb*sqrt(squeeze(sum(C.Vh(:, 2:end, :).^2, 2)));
else
  if disc
    load(sprintf('./DATS/LD_%s_NLINRESP_%s_%.2f_H%d.mat',nlin,mu,upordown,Nhmax), 'C', 'Ubws', 'h');
  else
    load(sprintf('./DATS/ROD_%s_NLINRESP_%s_%.2f_H%d.mat',nlin,mu,upordown,Nhmax), 'C', 'Ubws', 'h');
  end
end

if disc 
  save(sprintf('./DATS/LD_%s_%.2f_NLINRESP_%s_H%d.mat',nlin,mu,upordown,Nhmax), 'C', 'Ubws', 'h');
else 
  save(sprintf('./DATS/ROD_%s_%.2f_NLINRESP_%s_H%d.mat',nlin,mu,upordown,Nhmax), 'C', 'Ubws', 'h');
end

%% Plotting
figure(1)
% clf()
plot(C.Ws*60/2/pi, C.Udyn_amp(1, :), '-', 'LineWidth', 2); hold on
plot(C.Ws*60/2/pi, ones(size(C.Ws))*nl_els.pars(1), 'k--', 'LineWidth', 2);
## plot(C.Ws*60/2/pi, C.Udyn_amp(1:2, :), '-', 'LineWidth', 2);
set(gca, 'fontsize', 20)
xlabel('Spin Speed (RPM)')
ylabel('Displacement Amplitude (m)')

%% Depiction
if depict
  figure(2)
  clf()
  Wradpss = linspace(Wstart, Wend, 4);
  for iM=1:length(Wradpss)
    subplot(2, 2, iM)
    X = reshape(interp1(C.Ws, Ubws(1:end-1,:)', Wradpss(iM))', Nd, Nhc);
    DEPICTROTOR(Xs, Ris, Ros, Disc, Lb*X, h, 0);
    axis equal
    set(gca, 'fontsize', 20)
    title(sprintf('Spin Speed: %f RPM', Wradpss(iM)*60/2/pi))
  end
end

%% Lissajous
if lissajous
  figure(3)
  clf()
  Wradpss = linspace(C.Ws(10), C.Ws(end), 4);
  Nt = 128;
  t = linspace(0, 2*pi, Nt+1);
  for iM=1:length(Wradpss)
    subplot(2, 2, iM)
    X = reshape(interp1(C.Ws, Ubws(1:end-1,:)', Wradpss(iM))', Nd, Nhc);
    x = TIMESERIES_DERIV(Nt, h, (Rb*X)', 0);

    plot(x(:,1), x(:,2), 'o-'); hold on
    plot(cos(t)*nl_els.pars(1), sin(t)*nl_els.pars(1), 'k--', 'LineWidth', 2);
    
    axis equal
    set(gca, 'fontsize', 20)
    title(sprintf('Spin Speed: %f RPM', Wradpss(iM)*60/2/pi))
  end
end
