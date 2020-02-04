function [T,X,Xd,Xdd] = HHTA(M,C,K,FN,Fnl,X0,Xd0,t0,t1,dt,a,b,g)
%HHTA returns the HHT-Alpha time march
% USAGE:
%   [T,X,Xd] = HHTA(M,C,K,@(t) FN(t),X0,Xd0,t0,t1,h,a,g,b);
% INPUTS:
%   M       : NxN Intertia matrix
%   C       : NxN Proportional damping matrix
%   K       : NxN Stiffness matrix
%   FN      : Function handle for nonhomogeneous part returning 
%               Nx1 vector
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
    
    Z1 = M+(1+a)*g*dt*C+(1+a)*b*dt^2*K;
    Z2 = M-(1+a)*(1-g)*dt*C-(1+a)*(0.5-b)*dt^2*K;
    Z3 = (1+a)*dt*K;
    
    Z1 = inv(Z1);
    Xdd0 = M\(FN(t0)+Fnl(t0,X0,Xd0)-C*Xd0-K*X0);
    
    T = t0:dt:t1;       Nt = length(T);
    X = zeros(N,Nt);    Xd = X;             Xdd = X;
    X(:,1) = X0;        Xd(:,1) = Xd0;      Xdd(:,1) = Xdd0;
    
    for i=2:length(T)
        Xdd(:,i) = Z1*(Z2*Xdd(:,i-1) - Z3*Xd(:,i-1) + ...
            FN(T(i)+(1+a)*dt) - FN(T(i)) + ...
            Fnl(T(i)+(1+a)*dt,X(:,i-1),Xd(:,i-1)) - Fnl(T(i),X(:,i-1)+Xd(:,i-1)*dt,Xd(:,i-1)+Xdd(:,i-1)*dt) );
        Xd(:,i) = Xd(:,i-1) + dt*((1-g)*Xdd(:,i-1)+g*Xdd(:,i));
        X(:,i) = X(:,i-1) + dt*Xd(:,i-1) + dt^2*((0.5-b)*Xdd(:,i-1)+...
            b*Xdd(:,i));
    end
end