clc
clear all 
addpath('../ROUTINES')

disc = true;
analyze = false;
depict = true;
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
Wstart = 5*2*pi/60;  % 5 rpm
Wend = 250*2*pi/60;  % 250 rpm

Nd = size(Mb, 1);

ds = (Wend-Wstart)/100;
Nt = 128;
h = [0 1];
Nhc = sum(h==0) + 2*sum(h~=0);

Copt = struct('Nmax', 1000, 'Display', 1, 'angopt', 1e-6);
if analyze
  Ubws = CONTINUE(@(Uw) ROTOR_RESFUN(Uw, Mb, Cb, Gb, Kb, Fgb, FCb, FSb, h, Nt), zeros(Nd*Nhc, 1), ...
		  Wstart, Wend, ds, Copt);

  D1 = HARMONICSTIFFNESS(0, 1, 0, 1, h);
  C.Ws = Ubws(end, :);
  C.Uh = reshape(Ubws(1:end-1, :), Nd, Nhc, []);
  C.Vh = reshape(kron(D1, eye(Nd))*Ubws(1:end-1,:).*C.Ws, Nd, Nhc, []);
  C.Udyn_amp = Rb*sqrt(squeeze(sum(C.Uh(:, 2:end, :).^2, 2)));
  C.Vdyn_amp = Rb*sqrt(squeeze(sum(C.Vh(:, 2:end, :).^2, 2)));
else
  if disc
    load('./DATS/ROD_LINRESP.mat', 'C', 'Ubws');
  else
    load('./DATS/LD_LINRESP.mat', 'C', 'Ubws');
  end
end

if disc 
  save('./DATS/ROD_LINRESP.mat', 'C', 'Ubws');
else 
  save('./DATS/LD_LINRESP.mat', 'C', 'Ubws');
end

%% Plotting
figure(1)
clf()
plot(C.Ws*60/2/pi, C.Vdyn_amp(1, :), '-', 'LineWidth', 2);
## plot(C.Ws*60/2/pi, C.Udyn_amp(1:2, :), '-', 'LineWidth', 2);

%% Depiction
if depict
  figure(2)
  clf()
  [~, iM] = max(C.Vdyn_amp(1,:));

  DEPICTROTOR(Xs, Ris, Ros, Disc, Lb*C.Uh(:, :, iM), h);
  axis equal
end
