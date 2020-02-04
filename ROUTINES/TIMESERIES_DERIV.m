function [x_t, h] = TIMESERIES_DERIV(Nt,h,X0,ord)
% TIMESERIES_DERIV.m is a function that takes the fourier coefficients of a
% variable for specific harmonics and returns the variable in the time
% domain.
%
% USAGE:
%   x_t = TIMESERIES_DERIV(Nt,h,X0,ord)
% INPUTS:
%   Nt		: Number of time points (has to be even)
%   h   	: (Nhx1) a vector containing the harmonics being used to evaluate the
%           	  function in the time domain.
%   w   	: the frequency of which the system is being analyzed
%   X0  	: (2*Nh+1)xnd The fourier coefficients of the displacement of the Duffing Eqn
%   ord 	: Order of derivative
%  OUTPUTS:
%   x_t 	: (Ntxnd) the solution X0 in the time domain.
%   h		: List of harmonics (unchanged if 0 was not
%   		  there/was in the first position)    
%  EXAMPLE:
%  t = linspace(0, 2*pi, 129); t(end) = [];
%  y = sin(t) + 2*cos(2*t) + 3*sin(3*t);
%  h = [0; 1; 2; 3; 4];
%  yh = GETFOURIERCOEFF(h, y);
%  y_t = TIMESERIES_DERIV(128, h, yh, 0);
%  dydt_t = TIMESERIES_DERIV(128, h, yh, 1);
%  % ... so on
    % % IMPLEMENTATION WITH COMPULSORY REQUIREMENT OF ZERO HARMONIC
    % h = reshape(h,length(h),1);
    % nd = size(X0,2);
    % Nh = max(h);
        
    % X0full = zeros(2*Nh+1,nd);
    % if h(1)==0
    %     X0full(1,:) = X0(1,:);
    % end
    % X0full(2*h(2:end),:) = X0(2:2:end,:);
    % X0full(2*h(2:end)+1,:) = X0(3:2:end,:);
    % if Nt==0
    %     Nt = 2*Nh+2;
    % elseif Nt<2*Nh+1
    %     disp('Insufficient time points - truncating harmonics');
    %     Nht = floor(Nt/2-1);
    %     X0full = X0full(1:(2*Nht+1),:);
    %     Nh = Nht;
    %     Nt = 2*Nh+2;
    % else
    %     Nht = floor(Nt/2-1);
    %     X0full = [X0full(:,:);zeros(2*(Nht-Nh),nd)];
    %     Nh = Nht;
    %     Nt = 2*Nh+2;
    % end
    
    
    % D1 = zeros((2*Nh+1));
    % for k=1:Nh
    %     cosrows = 1 + (k-1)*2 + 1;
    %     sinrows = 1 + (k-1)*2 + 2;
    %     D1(cosrows,sinrows) = k;
    %     D1(sinrows,cosrows) = -k;
    % end
    % D = D1^ord;
    % X0full = D*X0full;
    
    
    % Xf = [2*X0full(1,:); X0full(2:2:end,:)-1j*X0full(3:2:end,:); zeros(1,nd); X0full((end-1):-2:2,:)+1j*X0full(end:-2:3,:);];
    % if size(Xf,1)~=Nt
    %     disp('error');
    % end
    % Xf = Xf*Nt/2;
    % x_t = real(ifft(Xf));
    
    % IMPLEMENTATION WITHOUT COMPULSORY REQUIREMENT OF ZERO
    % HARMONIC
    h = reshape(h,length(h),1);
    nd = size(X0,2);
    Nh = max(h);
    zi = find(h==0);
    h([1 zi]) = h([zi 1]);
    
    X0full = zeros(2*Nh+1,nd);
    if h(1)==0
        X0full(1,:) = X0(1,:);
        X0full(2*h(2:end),:) = X0(2:2:end,:);
        X0full(2*h(2:end)+1,:) = X0(3:2:end,:);
    else
        X0full(2*h(1:end),:) = X0(1:2:end,:);
        X0full(2*h(1:end)+1,:) = X0(2:2:end,:);
    end
    if Nt==0
        Nt = 2*Nh+2;
    elseif Nt<2*Nh+1
        disp('Insufficient time points - truncating harmonics');
        Nht = floor(Nt/2-1);
        X0full = X0full(1:(2*Nht+1),:);
        h(h>Nht) = [];
        Nh = Nht;
        Nt = 2*Nh+2;
    else
        Nht = floor(Nt/2-1);
        X0full = [X0full(:,:);zeros(2*(Nht-Nh),nd)];
        Nh = Nht;
        Nt = 2*Nh+2;
    end
    if ord>0
        D1 = zeros((2*Nh+1));
        for k=1:Nh
            cosrows = 1 + (k-1)*2 + 1;
            sinrows = 1 + (k-1)*2 + 2;
            D1(cosrows,sinrows) = k;
            D1(sinrows,cosrows) = -k;
        end
        D = D1^ord;
        X0full = D*X0full;
    end
    Xf = [2*X0full(1,:); X0full(2:2:end,:)-1j*X0full(3:2:end,:); zeros(1,nd); X0full((end-1):-2:2,:)+1j*X0full(end:-2:3,:);];
    if size(Xf,1)~=Nt
        disp('error');
    end
    Xf = Xf*Nt/2;
    x_t = real(ifft(Xf));    
end

%% VERIFIED
