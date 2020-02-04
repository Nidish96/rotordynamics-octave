function [] = DEPICTROTOR(Xs, Ris, Ros, Discs, varargin)
				%DEPICTROTOR depicts the rotor
  Nn = length(Xs);
  Ne = Nn-1;
  Nd = length(Discs);

  if nargin<7 || varargin{3}~=0
    colos = colormap(jet(Ne+Nd));
				% Rotating Elements
    for e=1:Ne
      [x, y, z] = cylinder([Ros(e) Ros(e)]);
      surf(z*(Xs(e+1)-Xs(e))+Xs(e), x, y, 'facecolor', colos(e,:)); hold on
      [r, t] = meshgrid([Ris(e) Ros(e)], linspace(0, 2*pi, 20));
      surf(Xs(e)*ones(size(r)), r.*cos(t), r.*sin(t), 'facecolor', colos(e,:));
      surf(Xs(e+1)*ones(size(r)), r.*cos(t), r.*sin(t), 'facecolor', colos(e,:));
      hold on
    end
				% Discs
    for d=1:length(Discs)
      [x, y, z] = cylinder([Discs(d).Rd Discs(d).Rd]);
      surf(z*Discs(d).t+Xs(Discs(d).nd)-Discs(d).t/2, x, y, 'facecolor', colos(Ne+d,:));
      [r, t] = meshgrid([Discs(1).Rs Discs(1).Rd], linspace(0, 2*pi, 20));
      surf((Xs(Discs(d).nd)-Discs(d).t/2)*ones(size(r)), r.*cos(t), r.*sin(t), 'facecolor', colos(Ne+d,:));
      surf((Xs(Discs(d).nd)+Discs(d).t/2)*ones(size(r)), r.*cos(t), r.*sin(t), 'facecolor', colos(Ne+d,:));
    end
  else
    plot3(Xs, Xs*0, Xs*0, 'ko-', 'LineWidth', 2); hold on
  end

  if length(varargin)==1
    Nt = 20;
    t = linspace(0, 2*pi, Nt);
    ux = varargin{1}(1:4:end, 1) + ...
         varargin{1}(1:4:end, 2)*cos(t) + ...
         varargin{1}(1:4:end, 3)*sin(t);
    tx = varargin{1}(3:4:end, 1) + ...
         varargin{1}(3:4:end, 2)*cos(t) + ...
         varargin{1}(3:4:end, 3)*sin(t);
    
    uy = varargin{1}(2:4:end, 1) + ...
         varargin{1}(2:4:end, 2)*cos(t) + ...
         varargin{1}(2:4:end, 3)*sin(t);
    ty = varargin{1}(4:4:end, 1) + ...
         varargin{1}(4:4:end, 2)*cos(t) + ...
         varargin{1}(4:4:end, 3)*sin(t);
    
    for n=1:Nn
      plot3(ones(Nt,1)*Xs(n), ux(n, :), uy(n, :), 'k-', ...
            'LineWidth', 2)
    end

    No = 10;

    Q = EBBM3D_ND2QP(diff(Xs), No);
    Xqp = Q(1:3:end, 1:5:end)*Xs(:);
    Q(1:3:end, :) = [];
    Q(:, 1:5:end) = [];
    Uqp = Q*reshape([ux(:,1)'; ty(:,1)'; uy(:,1)'; tx(:,1)'], [], 1);
    plot3(Xqp, Uqp(1:2:end), Uqp(2:2:end), 'k--', 'LineWidth', 2)
    plot3(Xs, ux(:, 1), uy(:, 1), 'ko', 'LineWidth', 2)
  elseif length(varargin)>=2  % Multiharmonic Input
    Nt = 128;
    t = linspace(0, 2*pi, Nt);
    h = varargin{2};
    ux = TIMESERIES_DERIV(Nt, h, varargin{1}(1:4:end, :)', 0)';
    tx = TIMESERIES_DERIV(Nt, h, varargin{1}(3:4:end, :)', 0)';

    uy = TIMESERIES_DERIV(Nt, h, varargin{1}(2:4:end, :)', 0)';
    ty = TIMESERIES_DERIV(Nt, h, varargin{1}(4:4:end, :)', 0)';

    for n=1:Nn
      plot3(ones(Nt,1)*Xs(n), ux(n, :), uy(n, :), 'k-', ...
            'LineWidth', 1)
    end

    No = 10;

    Q = EBBM3D_ND2QP(diff(Xs), No);
    Xqp = Q(1:3:end, 1:5:end)*Xs(:);
    Q(1:3:end, :) = [];
    Q(:, 1:5:end) = [];
    Uqp = Q*reshape([ux(:,1)'; ty(:,1)'; uy(:,1)'; tx(:,1)'], [], 1);
    plot3(Xqp, Uqp(1:2:end), Uqp(2:2:end), 'k--', 'LineWidth', 2)
    plot3(Xs, ux(:, 1), uy(:, 1), 'ko', 'LineWidth', 2)
  end
end
