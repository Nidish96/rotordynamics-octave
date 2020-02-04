function [T,X,Xd,Xdd] = HHTA_GEN(M,C,K,FN,nl_els,X0,Xd0,t0,t1,dt,a,b,g)
%HHTA returns the HHT-Alpha time march
% USAGE:
%   [T,X,Xd] = HHTA(M,C,K,@(t) FN(t),X0,Xd0,t0,t1,h,a,g,b);
% INPUTS:
%   M       : NxN Intertia matrix
%   C       : NxN Proportional damping matrix
%   K       : NxN Stiffness matrix
%   FN      : Function handle for nonhomogeneous part returning 
%               Nx1 vector
%   nl_els  : Array of structs with following options
% 	type 	 : supported: 'clearance', 'stator'
%	pars 	 : vector of parameters interpreted for each type as,
%			'clearance' -> [gap; stiffness]
%  			'stator' -> [gap; stiffness; mu]
%	sel_shape: NnlxNd matrix selecting non-linear DoFs for applying forces
%			'clearance' -> 2xNd (chosen DoFs: ux, uy)
%			'stator' -> 2xNd (chosen DoFs: ux, uy)
%	f_shape: NdxNnl matrix applying calculated forces to system
%			'clearance' -> Ndx2
%			'stator' -> Ndx2
%   X0      : Nx1 Vector of initial displacements
%   Xd0     : Nx1 Vector of initial velocities
%   t0      : Initial time
%   t1      : Final time
%   dt      : Time step
%   a       : Alpha parameter
%   b       : Beta parameter
%   g       : Gamma parameter
% OUTPUTS:
%   T       : 1xNt Time vector
%   X       : NxNt Displacement time series
%   Xd      : NxNt Velocity time series
%   Xdd     : NxNt Acceleration time series

  N = length(X0);
  X0 = reshape(X0,N,1);
  Xd0 = reshape(Xd0,N,1);
  Fnl_p = zeros(N, 1);  % "prev"
  Fnl_n = zeros(N, 1);  % "next"
  
  Z1 = M+(1+a)*g*dt*C+(1+a)*b*dt^2*K;
  Z2 = M-(1+a)*(1-g)*dt*C-(1+a)*(0.5-b)*dt^2*K;
  Z3 = (1+a)*dt*K;
  
  Z1 = inv(Z1);
  Fnl_p = NL_EVAL(t0, X0, Xd0, nl_els);
  Xdd0 = M\(FN(t0)-C*Xd0-K*X0 - Fnl_p);
## keyboard

  T = t0:dt:t1;       Nt = length(T);
  X = zeros(N,Nt);    Xd = X;             Xdd = X;
  X(:,1) = X0;        Xd(:,1) = Xd0;      Xdd(:,1) = Xdd0;
  
  for i=2:length(T)
    ## Fnl_n = NL_EVAL(T(i)+(1+a)*dt, X(:,i-1),              Xd(:,i-1),               nl_els);
    ## Fnl_p = NL_EVAL(T(i),          X(:,i-1)+Xd(:,i-1)*dt, Xd(:,i-1)+Xdd(:,i-1)*dt, nl_els);
    Fnl_n = NL_EVAL(T(i)+(1+a)*dt, X(:,i-1)+Xd(:,i-1)*dt*(2+a), Xd(:,i-1)+Xdd(:,i-1)*dt*(2+a), nl_els);
    Fnl_p = NL_EVAL(T(i),          X(:,i-1)+Xd(:,i-1)*dt,       Xd(:,i-1)+Xdd(:,i-1)*dt,       nl_els);
    
    Xdd(:,i) = Z1*(Z2*Xdd(:,i-1) - Z3*Xd(:,i-1) + ...
                   FN(T(i)+(1+a)*dt) - FN(T(i)) + ...
		   Fnl_p - Fnl_n );
    Xd(:,i) = Xd(:,i-1) + dt*((1-g)*Xdd(:,i-1)+g*Xdd(:,i));
    X(:,i) = X(:,i-1) + dt*Xd(:,i-1) + dt^2*((0.5-b)*Xdd(:,i-1)+...
                                             b*Xdd(:,i));

    fprintf('%d/%d %e\n', i, length(T), max(abs(X(:, i))));
    if ~isfinite(abs(X(:, i)))
      error('Out of Bounds!')
    end
  end
end

function [Fnl] = NL_EVAL(t, X, Xd, nl_els)
  Fnl = zeros(size(X));
  for iN=1:length(nl_els)
    q = nl_els(iN).sel_shape*X;
    qd = nl_els(iN).sel_shape*Xd;
    
    switch(nl_els(iN).type)
      case {'clearance', 'Clearance'}
	gap = nl_els(iN).pars(1);
	kr = nl_els(iN).pars(2);

	r = sqrt(sum(q.^2));
	if r~=0
	  fnl = kr*max(r-gap, 0).*(q./r);
	else
	  fnl = [0; 0];
	end
      case {'stator', 'Stator'}
	gap = nl_els(iN).pars(1);
	kr = nl_els(iN).pars(2);
	mu = nl_els(iN).pars(3);

	x = q(1); y = q(2); xd = qd(1); yd = qd(2);
	r = sqrt(sum(q.^2));
	fn = kr*max(r-gap, 0);

	if r~=0
	  fnl = fn.*(q./r) + mu*fn.*sign(x.*yd-xd.*y).*[-y;x]./r;
	else
	  fnl = [0; 0];
	end
    end
    ## if ~isempty(find(fnl~=0))
    ##   disp('Engaged!');
    ##   # keyboard
    ## end
    Fnl = Fnl + nl_els(iN).f_shape*fnl;
  end
end
