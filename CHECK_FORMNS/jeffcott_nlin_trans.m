#!/usr/bin/octave
clc
clear all
addpath('../ROUTINES');

nonlinearity = 'stator';
analyze = true;
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
nl_els = struct('type', nonlinearity, 'pars', [1e-6; 1e4; 0.5], ...
		'sel_shape', Lb((Disc.nd-1)*4+(1:2), :), ...
		'f_shape', Lb((Disc.nd-1)*4+(1:2), :)');

%% HHT-Alpha Integrator
ABG = [0, 1/4, 1/2];  % Unconditionally Stable Newmark-Alpha
## ABG = [0, 1/6, 1/2];  % Implicit linear acceleration
## ABG = [-0.1, 1/6, 1/2];  % HHT-Alpha

iopts = struct('a', ABG(1), 'b', ABG(2), 'g', ABG(3), 'ITMAX', 100, ...
	       'ITOPT', 10, 'etol', 1e-3, 'reletol', 1e-6, 'rtol', 1e-6, 'utol', 1e-6, ...
	       'Display', false, 'adapt', false);

Tmax = 10;
Tcut = 8;
fs = 2^14;

Nd = size(Mb, 1);

X0 = zeros(Nd, 1);
Xd0 = zeros(Nd, 1);

Wrpms = linspace(100, 2e5, 100);

## Wrpms = 9e4;
fid = fopen('log.dat', 'w+');
fprintf(fid, 'Log file\n');
if analyze
  fprintf('+----------------- Starting Transient Analysis ----------------+\n');
  for i = 1:100
    Wrpm = Wrpms(i);
    Wradps = Wrpm*2*pi/60;

    Fn = @(t) Wradps^2*(FCb.*cos(Wradps*t)+FSb.*sin(Wradps*t))+Fgb;

    try
      tic
      [T, X, Xd, Xdd] = HHTA_NL(Mb, Cb-Wradps*Gb, Kb, Fn, nl_els, X0, Xd0, 0, Tmax, 1/fs, iopts);
      toc
      save(sprintf('./DATS/jeffcott_%s_nlin_trans_%d.mat',nonlinearity, i), 'T', 'X', 'Xd', 'Xdd', 'Wrpm', 'Wradps');
    catch me
      fprintf(fid, '%d %frpm\n', i, Wrpm);
      disp('not complete');
    end

    fprintf('Done %d/%d\n', i, length(Wrpms))
  end
end
fclose(fid);

% Single Data Plotting
##i = 1;
##load(sprintf('./DATS/jeffcott_%s_nlin_trans_%d.mat', nonlinearity, i));
        
figure(1)
clf()
plot(T, Rb(1:2,:)*X, '.-', 'LineWidth', 2)
grid on
xlabel('Time (s)')
ylabel('Transverse Deflections (m)')
legend('Horizontal', 'Vertical')


ti = find((T(1:end-1)-Tcut).*(T(2:end)-Tcut)<=0); ti = min(ti);
[freq, Xf] = FFTFUN(T(ti:end)-T(ti), (Rb(1:2, :)*X(:, ti:end))');

figure(2)
clf()
subplot(2,1,1)
semilogy(freq(1:size(Xf,1)), abs(Xf), 'LineWidth', 2)
grid on
xlabel('Frequency (Hz)')
ylabel('Transverse Deflections (m)')
legend('Horizontal', 'Vertical')
xlim([0 fs/2])

subplot(2,1,2)
plot(freq(1:size(Xf,1)), rad2deg(angle(Xf)), 'LineWidth', 2)
grid on
xlabel('Frequency (Hz)')
ylabel('Phase (degs)')
xlim([0 fs/2])

figure(3)

No = 10;
Q = EBBM3D_ND2QP(diff(Xs), No);
Xqp = Q(1:3:end, 1:5:end)*Xs(:);
Q(1:3:end, :) = [];
Q(:, 1:5:end) = [];

Nt = length(T);
Np = 1000;

xm = max(max(abs(X(1:4:end,:))));
ym = max(max(abs(X(2:4:end,:))));
for i=fix(linspace(1, Nt, Np))
  clf()
  plot3(Xs, Xs*0, Xs*0); hold on
  plot3(Xs(Disc.nd)*ones(Nt, 1), nl_els.pars(1)*cos(T), nl_els.pars(1)*sin(T), 'k--')

  uqp = Q*Lb(reshape([(1:4:Nn*4)' (3:4:Nn*4)' (2:4:Nn*4)' (4:4:Nn*4)']', [], 1), :)*X(:, i);
  
  plot3(Xqp, uqp(1:2:end), uqp(2:2:end), 'k-', 'LineWidth', 2);
  plot3(Xs, Lb(1:4:end, :)*X(:, i), Lb(2:4:end, :)*X(:, i), 'ko', 'LineWidth', 2);

  grid on;
  xlim(Xs([1 end]))
  ylim(xm*[-1 1])
  zlim(ym*[-1 1])
  title(sprintf('Spin Speed: %.2f rpm\t Frame %d/%d', Wradps*60/2/pi, i, Nt))
  pause(0.0001)
end
# legend('Horizontal', 'Vertical')

[S, f, t] = specgram(Rb(1,:)*X, 256, fs);
figure(4); imagesc(t, f, log(abs(S))); set (gca, "ydir", "normal");

##				% Multiple Data Plotting
##% Spectrogram Evolution
##for i=fix(linspace(1, length(Wrpms), 10))
##  load(sprintf('./DATS/jeffcott_clearance_nlin_trans_%d.mat',i))
##  
##  figure(1)
##  clf()
##  
##  specgram(Rb(1,:)*X, 256, fs);
##  title(sprintf('Spin Speed = %f', Wrpms(i)))
##
##  pause(0.0001)
##  fprintf('Plotted %d\n', i);
##end
##
##				% Waterfall plot
##figure(3)
##clf()
##for i=1:length(Wrpms)
##  load(sprintf('./DATS/jeffcott_clearance_nlin_trans_%d.mat',i));
##  ti = find((T(1:end-1)-Tcut).*(T(2:end)-Tcut)<=0); ti = min(ti);
##  [freq, Xf] = FFTFUN(T(ti:end)-T(ti), (Rb(1:2, :)*X(:, ti:end))');
##
##  plot3(freq, Wrpms(i)*ones(size(freq)), abs(Xf(:, 1)), 'k-', 'LineWidth', 2); hold on
##  pause(0.0001)
##end
##xlabel('Frequency (Hz)')
##ylabel('Spin Speed (rpm)')
##zlabel('Amplitude (m)')
