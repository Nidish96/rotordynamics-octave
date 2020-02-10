function [T,X,Xd,Xdd] = HHTA_NL(M,C,K,FN,nl_els,X0,Xd0,t0,t1,dt,opts)
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

  %% Interpret Options
  a = opts.a;
  b = opts.b;
  g = opts.g;
  
  N = length(X0);
  X0 = reshape(X0,N,1);
  Xd0 = reshape(Xd0,N,1);
  Fnl_p = zeros(N, 1);  % "prev"
  Fnl_n = zeros(N, 1);  % "next"
  Jnl_n = zeros(N);
  Fl_p = zeros(N, 1);
  Fl_n = zeros(N, 1);
  R = zeros(N, 1);
  J = zeros(N);
  
  Z1 = M+(1+a)*g*dt*C+(1+a)*b*dt^2*K;
  Z2 = M-(1+a)*(1-g)*dt*C-(1+a)*(0.5-b)*dt^2*K;
  Z3 = (1+a)*dt*K;
  
  Fnl_p = NL_EVAL(t0, X0, Xd0, nl_els);
  Xdd0 = M\(FN(t0)-C*Xd0-K*X0 - Fnl_p);
  
  Ntmax = fix((t1-t0)/dt+10);
  dt0 = dt;
  
  T = zeros(Ntmax, 1);  T(1) = t0;
  X = zeros(N,Ntmax);    Xd = X;             Xdd = X;
  X(:,1) = X0;        Xd(:,1) = Xd0;      Xdd(:,1) = Xdd0;

  i = 2;
  flag = 0;
  while T(i-1)<t1 && i<=Ntmax
    T(i) = T(i-1) + dt;
				% Explicit Predictor: Constant acceleration 
    Fl_p = FN(T(i-1));
    Fnl_p = NL_EVAL(T(i-1), X(:,i-1), Xd(:, i-1), nl_els);
    Xdd(:, i) = Xdd(:, i-1);  % Constant Acceleration
    
				% Corrector Steps
    Fl_n = FN(T(i-1)+dt*(1+a));
    [Fnl_n, Jnl_n] = NL_EVAL(T(i-1)+dt*(1+a), ...
			     X(:, i-1)+(1+a)*dt*Xd(:,i-1)+(1+a)*dt^2*((0.5-b)*Xdd(:, i-1)+b*Xdd(:, i)), ...
			     Xd(:, i-1)+(1+a)*dt*((1-g)*Xdd(:, i-1)+g*Xdd(:, i)), ...
			     nl_els);
    R = Z1*Xdd(:, i) - Z2*Xdd(:, i-1) + Z3*Xd(:, i-1) + Fnl_n - Fnl_p - Fl_n + Fl_p;
    J = Z1 + Jnl_n*(1+a)*dt^2*b;
    du = -J\R;
    e0 = abs(R'*du);
    r0 = mean(R.^2);
    u0 = mean(du.^2);
    
    e = e0;
    r = r0;
    u = u0;
    it = 0;

    flag = 1*(e/e0<opts.reletol) + 2*(e<opts.etol) + 4*(r<opts.rtol) + 8*(u<opts.utol);
    if opts.Display
      fprintf('ITN, E, E/E0, r, du\n%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, flag);
    end

    flag = 0;
    while (flag < 3)
      Xdd(:, i) = Xdd(:, i) + du;
      it = it+1;

      [Fnl_n, Jnl_n] = NL_EVAL(T(i-1)+dt*(1+a), ...
			       X(:, i-1)+(1+a)*dt*Xd(:,i-1)+(1+a)*dt^2*((0.5-b)*Xdd(:, i-1)+b*Xdd(:, i)), ...
			       Xd(:, i-1)+(1+a)*dt*((1-g)*Xdd(:, i-1)+g*Xdd(:, i)), ...
			       nl_els);
      R = Z1*Xdd(:, i) - Z2*Xdd(:, i-1) + Z3*Xd(:, i-1) + Fnl_n - Fnl_p - Fl_n + Fl_p;
      J = Z1 + Jnl_n*(1+a)*dt^2*b;
      du = -J\R;
      e = abs(R'*du);
      r = mean(R.^2);
      u = mean(du.^2);
      
      flag = 1*(e/e0<opts.reletol) + 2*(e<opts.etol) + 4*(r<opts.rtol) + 8*(u<opts.utol);
      if opts.Display
	fprintf('%d, %e, %e, %e, %e: %d\n', it, e, e/e0, mean(R.^2), mean(du.^2), flag);
      end
      
      if it>opts.ITMAX
	if flag == 0
	  disp('No convergence: trying to halve size if allowed');

	  if opts.adapt
	    if dt>dt0/5
	      break
	    else
	      error('nah');
	    end
	  else
	    error('nah');
	  end
	else
	  flag = 3;
	end
      end      
    end
    if opts.Display
      fprintf('\n');
    end
    if opts.adapt
      if it>opts.ITOPT
	dt = dt/2;
      elseif it<opts.ITOPT && dt<dt0*2
	dt = dt*2;
      end
    end
    
    Xd(:,i) = Xd(:,i-1) + dt*((1-g)*Xdd(:,i-1)+g*Xdd(:,i));
    X(:,i) = X(:,i-1) + dt*Xd(:,i-1) + dt^2*((0.5-b)*Xdd(:,i-1)+b*Xdd(:,i));

    if opts.Display
      fprintf('%.4e/%.4e %.4e\n', T(i), t1, dt);
    end
				% Step Successful
    i = i+1;
    if ~isfinite(abs(X(:, i)))
      error('Out of Bounds!')
    end
  end

  T = T(1:i-1);
  X = X(:, 1:i-1);  Xd = Xd(:, 1:i-1);  Xdd = Xdd(:, 1:i-1);
end

function [Fnl, Jnl] = NL_EVAL(t, X, Xd, nl_els)
  Fnl = zeros(size(X));
  Jnl = zeros(length(X));
  for iN=1:length(nl_els)
    q = nl_els(iN).sel_shape*X;
    qd = nl_els(iN).sel_shape*Xd;
    
    switch(nl_els(iN).type)
      case {'clearance', 'Clearance'}
        gap = nl_els(iN).pars(1);
        kr = nl_els(iN).pars(2);

        x = q(1); y = q(2); xd = qd(1); yd = qd(2);
        r = sqrt(sum(q.^2));
        if r~=0
          fnl = kr*max(r-gap, 0).*(q./r);

          jnl_xx = kr*((r-gap)./r+gap*x.^2./r.^3);
          jnl_xy = kr*gap.*x.*y./r.^3;

          jnl_yx = kr*gap.*x.*y./r.^3;
          jnl_yy = kr*((r-gap)./r+gap*y.^2./r.^3);
        else
          fnl = [0; 0];

          jnl_xx = 0;
          jnl_xy = 0;

          jnl_yx = 0;
          jnl_yy = 0;
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

          jnl_xx = kr*((r-gap)./r+gap*x.^2./r.^3) + ...
             mu*(kr*((r-gap)./r+gap*x.^2./r.^3)).*sign(x.*yd-xd.*y).*(-y./r) + ...
             mu*fn.*sign(x.*yd-xd.*y).*(y.*x./r.^3);
          jnl_xy = (kr*gap.*x.*y./r.^3) + ...
             mu*((kr*gap.*x.*y./r.^3)).*sign(x.*yd-xd.*y).*(-y./r) + ...
             mu*fn.*sign(x.*yd-xd.*y).*(-x.^2./r.^3);

          jnl_yx = (kr*gap.*x.*y./r.^3) + ...
             mu*((kr*gap.*x.*y./r.^3)).*sign(x.*yd-xd.*y).*x./r + ...
             mu*fn.*sign(x.*yd-xd.*y).*(y.^2./r.^3);
          jnl_yy = kr*((r-gap)./r+gap*y.^2./r.^3) + ...
             mu*(kr*((r-gap)./r+gap*y.^2./r.^3)).*sign(x.*yd-xd.*y).*(x./r) + ...
             mu*fn.*sign(x.*yd-xd.*y).*(-x.*y./r.^3);
        else
          fnl = [0; 0];

          jnl_xx = 0;
          jnl_xy = 0;

          jnl_yx = 0;
          jnl_yy = 0;	  
        end
      case {'Cubic', 'cubic'}
        bx = nl_els(iN).pars(1);
        by = nl_els(iN).pars(2);
        fnl = [bx;by].*q.^3;

        jnl_xx = 3*bx*q(1)^2;
        jnl_xy = 0;

        jnl_yx = 0;
        jnl_yy = 3*by*q(2)^2;
      case {'Test', 'test'}
        kx = nl_els(iN).pars(1);
        ky = nl_els(iN).pars(2);
        fnl = [kx;ky].*q;

        jnl_xx = kx;
        jnl_xy = 0;

        jnl_yx = 0;
        jnl_yy = ky;
      otherwise
	      error ('Unimplemented');
    end
    Fnl = Fnl + nl_els(iN).f_shape*fnl;
    Jnl = Jnl + nl_els(iN).f_shape(:, 1)*jnl_xx*nl_els(iN).sel_shape(1,:) + ...
	  nl_els(iN).f_shape(:, 1)*jnl_xy*nl_els(iN).sel_shape(2,:) + ...
	  nl_els(iN).f_shape(:, 2)*jnl_yx*nl_els(iN).sel_shape(1,:) + ...
	  nl_els(iN).f_shape(:, 2)*jnl_yy*nl_els(iN).sel_shape(2,:);
  end
end
