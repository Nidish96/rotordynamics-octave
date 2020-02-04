clc
clear all
addpath('../ROUTINES')

%% Parameters
E = 210e9;
rho = 7800;

Xs  = [0, 10, 20, 25, 27.5, 30]*1e-3;
Ris = [0,  0,  0,  0,   0 ,  0]; 
Ros = [1,  1,  1,  1,   1 ,  1]*1e-3;

Ros_wdisc = [1,  1,  1,  1,   4 ,  4]*1e-3;

Nn = length(Xs);
Ne = Nn-1;

g = 9.81;

%% System without Disc in FE model
WOD.M = zeros(Nn*4);
WOD.K = zeros(Nn*4);
WOD.G = zeros(Nn*4);

Me = zeros(8);
Ge = zeros(8);
Ke = zeros(8);

for e=1:Ne
    [Me, Ge, Ke, ~] = ROTELMATS(E, rho, Xs(e+1)-Xs(e), Ris(e), Ros(e), g);
    
    is = (e-1)*4+(1:8);
    WOD.M(is, is) = WOD.M(is, is) + Me;
    WOD.G(is, is) = WOD.G(is, is) + Ge;
    WOD.K(is, is) = WOD.K(is, is) + Ke;
end

% Discrete Models of Disc
Disc = DISC(Nn-1, Ros(end), 4e-3, Xs(end)-Xs(end-2), rho, 0, 0);
[Md, Gd] = Disc.MATS();
WOD.M((Disc.nd-1)*4+(1:4), (Disc.nd-1)*4+(1:4)) = WOD.M((Disc.nd- ...
                                                  1)*4+(1:4), ...
                                                  (Disc.nd-1)*4+(1:4)) ...
    + Md;
WOD.G((Disc.nd-1)*4+(1:4), (Disc.nd-1)*4+(1:4)) = WOD.G((Disc.nd- ...
                                                  1)*4+(1:4), ...
                                                  (Disc.nd-1)*4+(1:4)) ...
    + Gd;

% Boundary Conditions
WOD.Lb = eye(Nn*4); WOD.Lb(:, 1:4) = [];
WOD.Mb = WOD.Lb'*WOD.M*WOD.Lb;
WOD.Gb = WOD.Lb'*WOD.G*WOD.Lb;
WOD.Kb = WOD.Lb'*WOD.K*WOD.Lb;

%% System with Modeled Disc
WD.M = zeros(Nn*4);
WD.K = zeros(Nn*4);
WD.G = zeros(Nn*4);

for e=1:Ne
  [Me, Ge, Ke, ~] = ROTELMATS(E, rho, Xs(e+1)-Xs(e), Ris(e), Ros_wdisc(e), g);

  is = (e-1)*4+(1:8);
  WD.M(is, is) = WD.M(is, is) + Me;
  WD.G(is, is) = WD.G(is, is) + Ge;
  WD.K(is, is) = WD.K(is, is) + Ke;
end

% Boundary Conditions
WD.Lb = eye(Nn*4); WD.Lb(:, 1:4) = [];
WD.Mb = WD.Lb'*WD.M*WD.Lb;
WD.Gb = WD.Lb'*WD.G*WD.Lb;
WD.Kb = WD.Lb'*WD.K*WD.Lb;


%% QEPS
WOD.Z = polyeig(WOD.Kb, -WOD.Gb, WOD.Mb);
[~, si] = sort(abs(WOD.Z));  WOD.Z = WOD.Z(si);

WD.Z = polyeig(WD.Kb, -WD.Gb, WD.Mb);
[~, si] = sort(abs(WD.Z));  WD.Z = WD.Z(si);

abs([imag([WOD.Z(1:2:end) WD.Z(1:2:end)]) WOD.Z(1:2:end)./WD.Z(1:2:end)])
