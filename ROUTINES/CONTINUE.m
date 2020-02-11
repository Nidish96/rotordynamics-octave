function [U, R, J] = CONTINUE(func, u0, lam0, lam1, ds, varargin)
%CONTINUE Conducts the continuation and solves the system
%
% USAGE:
% INPUTS:
% OUTPUTS:

				% Default options
  Copt = struct('Nmax', 100, 'tole', 1e-6, 'tolr', 1e-6, 'ITMAX', 1000, ...
                'dsmax', ds*5, 'dsmin', ds/5, 'angopt', pi/6, 'startdir', 1,...
                'Display', 1, 'nev', 1, 'adj', 1,...
                'opt', optimset('fsolve'));
  if nargin==6
    nflds = fieldnames(varargin{1});
    for i=1:length(nflds)
      Copt.(nflds{i}) = varargin{1}.(nflds{i});
    end
  end

				% Allocations
  U = zeros(length(u0)+1, Copt.Nmax);
  R = zeros(length(u0), Copt.Nmax);
  J = zeros(length(u0), length(u0)+1, Copt.Nmax);
  Jtmp = zeros(length(u0)+1);
  Rtmp = zeros(length(u0)+1, 1);
  
				% Correct initial solution
  [u0s, R(:, 1), eflag] = fsolve(@(u) func([u; lam0]), u0, Copt.opt);
 %     [u0s, R(:, 1), eflag] = NSOLVE(@(u) func([u; lam0]), u0, Copt);
  if eflag < 0
    error('Initial point non-convergent!');
  elseif Copt.Display
    disp('Initial Point Converged')
  end
  U(:, 1) = [u0s; lam0];
  [~, dRdX0, dRdlam0] = func(U(:, 1));
  J(:, :, 1) = [dRdX0 dRdlam0];
  
				% BEGIN CONTINUATION
  lam  = lam0;  % current parameter
  lamp = lam0;  % previous parameter
  n    = 1;
  
				% Initial tangent
  z = -dRdX0\dRdlam0;
  z(~isfinite(z)) = 0;
  al = 1.0/sqrt(1+z'*z)*Copt.startdir;
  
  alp = al;
  zp  = z;
  
  duds   = al*[z; 1];
  uguess = U(:, 1) + ds*duds;
  while ( (lam-lam1)*(lamp-lam1) >= 0 && n<Copt.Nmax )
    [U(:, n+1), Rtmp, eflag, out, Jtmp] = fsolve(@(u) EXRES(func, u, U(:, n), duds, ds), uguess, Copt.opt);
%         [U(:, n+1), Rtmp, eflag, out, Jtmp] = NSOLVE(@(u) EXRES(func, u, U(:, n), duds, ds), uguess, Copt);
    if eflag<=0
      if ds == Copt.dsmin
        if max(abs(uguess))<eps
          disp('Diverged!');
          break;
        else
          disp('Diverged! Trying with first initial guess!');
          uguess = U(:,1);
          continue;
        end
      else
        disp('Diverged - reducing step-size');
        ds = Copt.dsmin;
        uguess = U(:, n) + ds*duds;
        continue;
      end
    end
    R(:, n+1)    = Rtmp(1:end-1);
    J(:, :, n+1) = Jtmp(1:end-1, :);
    
				% Step
    lamp = lam;
    lam  = U(end, n+1);
    
				% tangent and predictor
    z = -J(:, 1:end-1, n+1)\J(:, end, n+1);
    z(~isfinite(z)) = 0;
    al = 1.0/sqrt(1+z'*z) * sign(alp*(1+zp'*z));
    
				% step size adaptation
    theta = acos(alp*al*(1+zp'*z));
    if theta>Copt.angopt % optimal angle
      ds = max(Copt.dsmin, ds/2);
    elseif out.iterations<=10
      ds = min(Copt.dsmax, ds*2);
    end
    if Copt.Display
      fprintf('%d %f %f %e %d\n', n+1, U(end,n+1), ds, theta, ...
              eflag);
    end
    
    alp = al;
    zp  = z;
    
    n = n+1;
    duds   = al*[z; 1];
    uguess = U(:, n) + ds*duds;
  end
  U = U(:, 1:n);
  R = R(:, 1:n);
  J = J(:, :, 1:n);
  
  if (lam-lam1)*(lamp-lam1) < 0
      disp('Continuation completed successfully');
  else 
      disp('Premature Termination');
  end
end

function [Re, Je] = EXRES(func, u, u0, up0, ds)
%EXRES is the residual function along with the continuation constraint
  [Re, dRdUe, dRdLe] = func(u);
  Re = [Re; up0'*(u-u0)-ds];
  Je = [dRdUe dRdLe; up0'];
end
