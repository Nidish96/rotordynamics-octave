function [R, dRdU, dRdw] = ROTOR_RESFUN(Uw, M, C, G, K, F0, Fc, Fs, h, ...
					Nt, varargin)
%ROTOR_RESFUN returns the residue function and its derivatives
%  USAGE:
%    [R, dRdU, dRdw] = ROTOR_RESFUN(Uw, M, C, G, K, L, F0, Fc, Fs, h, Nt, nl_els)
%  INPUTS:
%    Uw 	: (Nd*Nhc+1)x1 vector of unknowns
%    M, C, G, K : NdxNd Rotordynamic matrices
%    F0         : Ndx1 zero-harmonic "static" forces
%    Fc, Fs     : Ndx1 centrifugal harmonic forces (cosine & sine) - a scaling of w^2 will be applied
%    h 		: List of harmonics (e.g. [0 1 2 3])
%    Nt		: AFT time steps
%    nl_els     : Array of structs with following options
% 		type 	 : supported: 'clearance', 'stator'
%		pars 	 : vector of parameters interpreted for each type as,
%				'clearance' -> [gap; stiffness]
%  				'stator' -> [gap; stiffness; mu]
%		sel_shape: NnlxNd matrix selecting non-linear DoFs for applying forces
%				'clearance' -> 2xNd (chosen DoFs: ux, uy)
%				'stator' -> 2xNd (chosen DoFs: ux, uy)
%		f_shape: NdxNnl matrix applying calculated forces to system
%				'clearance' -> Ndx2
%				'stator' -> Ndx2
%  OUTPUTS:
%    R		: (Nd*Nhc)x1 residue vector
%    dRdU	: (Nd*Nhc)x(Nd*Nhc) residue jacobian
%    dRdw 	: (Nd*Nhc)x1 spin speed derivative

  %% Lower Case: time domain; Upper case: frequency domain
  
  if length(varargin)==1
    nl_els = varargin{1};
  else
    nl_els = [];
  end

  Nd = size(M,1);
  Nhc = sum(h==0)+2*sum(h~=0);
  
  w    = Uw(end);
  E    = HARMONICSTIFFNESS(M, C-w*G, K, w, h);
  dEdw = HARMONICSTIFFNESS(2*M/w, C/2-2*G, zeros(size(K)), w, h);
  D1   = HARMONICSTIFFNESS(0, 1, 0, w, h);

				% Non-Linearities
  t = linspace(0, 2*pi, Nt+1); t(end) = [];
  cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);  % Relevant Harmonic basis functions in columns
  sct = TIMESERIES_DERIV(Nt, h, D1, 0);  % Basis derivatives
  
  FNL = zeros(Nd*Nhc, 1);
  JNL = zeros(Nd*Nhc); 
  for iN=1:length(nl_els)
    Q = (nl_els(iN).sel_shape*reshape(Uw(1:Nd*Nhc), Nd, Nhc))';  % Harmonics of non-linear DoFs in columns
    q = TIMESERIES_DERIV(Nt, h, Q, 0);  % Time series of non-linear DoFs in columns
    qdot = TIMESERIES_DERIV(Nt, h, Q, 1);  % Time series of non-linear DoFs velocities in columns

    switch(nl_els(iN).type)
      case {'clearance', 'Clearance'}
	gap = nl_els(iN).pars(1);
	kr = nl_els(iN).pars(2);

	x = q(:, 1);
	y = q(:, 2);
	r = sqrt(sum(q.^2, 2));

	fnl = kr*max(r-gap, 0).*(q./r);
	jnl_xx = kr*((r-gap)./r+gap*x.^2./r.^3).*(r>gap);
	jnl_xy = (kr*gap.*x.*y./r.^3).*(r>gap);
	
	jnl_yx = (kr*gap.*x.*y./r.^3).*(r>gap);
	jnl_yy = kr*((r-gap)./r+gap*y.^2./r.^3).*(r>gap);

	%% Avoid divisions by r=0
	fnl(r==0, :) = 0;
	jnl_xx(r==0, :) = 0;
	jnl_xy(r==0, :) = 0;
	jnl_yx(r==0, :) = 0;
	jnl_yy(r==0, :) = 0;
      case {'stator', 'Stator'}
	gap = nl_els(iN).pars(1);
	kr = nl_els(iN).pars(2);
	mu = nl_els(iN).pars(3);

	x = q(:, 1);
	y = q(:, 2);
	r = sqrt(sum(q.^2, 2));
	xd = qdot(:, 1);
	yd = qdot(:, 2);

	fn = kr*max(r-gap, 0);

	fnl = fn.*(q./r) + ...
	      mu*fn.*sign(x.*yd-xd.*y).*[-y x]./r;

	jnl_xx = kr*((r-gap)./r+gap*x.^2./r.^3).*(r>gap) + ...
		 mu*(kr*((r-gap)./r+gap*x.^2./r.^3).*(r>gap)).*sign(x.*yd-xd.*y).*(-y./r) + ...
		 mu*fn.*sign(x.*yd-xd.*y).*(y.*x./r.^3);
	jnl_xy = (kr*gap.*x.*y./r.^3).*(r>gap) + ...
		 mu*((kr*gap.*x.*y./r.^3).*(r>gap)).*sign(x.*yd-xd.*y).*(-y./r) + ...
		 mu*fn.*sign(x.*yd-xd.*y).*(-x.^2./r.^3);

	jnl_yx = (kr*gap.*x.*y./r.^3).*(r>gap) + ...
		 mu*((kr*gap.*x.*y./r.^3).*(r>gap)).*sign(x.*yd-xd.*y).*x./r + ...
		 mu*fn.*sign(x.*yd-xd.*y).*(y.^2./r.^3);
	jnl_yy = kr*((r-gap)./r+gap*y.^2./r.^3).*(r>gap) + ...
		 mu*(kr*((r-gap)./r+gap*y.^2./r.^3).*(r>gap)).*sign(x.*yd-xd.*y).*(x./r) + ...
		 mu*fn.*sign(x.*yd-xd.*y).*(-x.*y./r.^3);

	%% Avoid divisions by r=0
	fnl(r==0, :) = 0;
	jnl_xx(r==0, :) = 0;
	jnl_xy(r==0, :) = 0;
	jnl_yx(r==0, :) = 0;
	jnl_yy(r==0, :) = 0;
    end
    Fnl = GETFOURIERCOEFF(h, fnl);
    Jnl_XX = GETFOURIERCOEFF(h, jnl_xx.*cst);
    Jnl_XY = GETFOURIERCOEFF(h, jnl_xy.*cst);

    Jnl_YX = GETFOURIERCOEFF(h, jnl_yx.*cst);
    Jnl_YY = GETFOURIERCOEFF(h, jnl_yy.*cst);

    FNL = FNL + reshape(nl_els(iN).f_shape*Fnl', Nd*Nhc, 1);

    JNL = JNL + kron(Jnl_XX, nl_els(iN).f_shape(:, 1).*nl_els(iN).sel_shape(1, :)) + ...
	  kron(Jnl_XY, nl_els(iN).f_shape(:, 1).*nl_els(iN).sel_shape(2, :)) + ...
	  kron(Jnl_YX, nl_els(iN).f_shape(:, 2).*nl_els(iN).sel_shape(1, :)) + ...
	  kron(Jnl_YY, nl_els(iN).f_shape(:, 2).*nl_els(iN).sel_shape(2, :));
  end
				% Residual
  R = E*Uw(1:Nd*Nhc) + FNL;
  R(1:Nd) = R(1:Nd) - F0;
  R(Nd+(1:2*Nd)) = R(Nd+(1:2*Nd)) - w^2*[Fc; Fs];

				% U-Jacobian
  dRdU = E + JNL;

				% w-Jacobian
  dRdw = dEdw*Uw(1:Nd*Nhc);
  dRdw(Nd+(1:2*Nd)) = dRdw(Nd+(1:2*Nd)) - 2*w*[Fc; Fs];
end
