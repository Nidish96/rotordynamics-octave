function [Me,Ge,Ke,Fge,MTe,MRe] = ROTELMATS(E, rho, L, Ri, Ro, g)
%ROTELMATS Returns the element matrices of the rotor element based on
%Euler-Bernouilli Beam Theory
% 
% USAGE: 
%   [Me,Ge,Ke,MTe,MRe] = ROTELMATS(E, rho, L, Ri, Ro);
% INPUTS:
%   E   : Youngs modulus
%   rho : Density
%   L   : Element Length
%   Ri  : Inner radius
%   Ro  : Outer Radius
% OUTPUTS:
%   Me  :
%   Ge  :
%   Ke  :
%   MTe :
%   MRe :

    A  = pi*(Ro^2-Ri^2);   % Area 
    m  = rho*A;            % Mass per unit length
    I  = pi*(Ro^4-Ri^4)/4; % Bending moment
    j  = I*rho;            % Bending M.I.
    jz = j*2;              % Polar M.I

    MTe = zeros(8);
    MRe = zeros(8);
    Ge  = zeros(8);
    Ke  = zeros(8);
    
    MTe(1,1:end) = [156  ,    0 ,  0,  22*L,   54  ,   0 ,  0, -13*L];
    MTe(2,2:end) = [156  , -22*L,  0,    0 ,   54  , 13*L,  0];
    MTe(3,3:end) = [4*L^2,    0 ,  0, -13*L, -3*L^2,   0 ];
    MTe(4,4:end) = [4*L^2,  13*L,  0,    0 , -3*L^2];
    MTe(5,5:end) = [156  ,    0 ,  0, -22*L];
    MTe(6,6:end) = [156  ,  22*L,  0];
    MTe(7,7:end) = [4*L^2,    0 ];
    MTe(8,8:end) = [4*L^2];
    MTe = (MTe + MTe') - diag(diag(MTe));
    MTe = MTe*m*L/420;  % Checked
    
    MRe(1,1:end) = [36   ,   0 , 0, 3*L, -36 ,   0 , 0, 3*L];
    MRe(2,2:end) = [36   , -3*L, 0,  0 , -36 , -3*L, 0];
    MRe(3,3:end) = [4*L^2,   0 , 0, 3*L, -L^2,   0 ];
    MRe(4,4:end) = [4*L^2, -3*L, 0,  0 , -L^2];
    MRe(5,5:end) = [36   ,   0 , 0, -3*L];
    MRe(6,6:end) = [36   ,  3*L, 0];
    MRe(7,7:end) = [4*L^2, 0];
    MRe(8,8:end) = [4*L^2]; 
    MRe = (MRe + MRe') - diag(diag(MRe));
    MRe = MRe*j/(30*L); % Checked
    
    Me = MTe+MRe;

    Ge(1,2:end) = [-36   ,  3*L,   0 , 0,  36, 3*L, 0];
    Ge(2,3:end) = [0     ,  3*L, -36 , 0,  0 , 3*L];
    Ge(3,4:end) = [-4*L^2,  3*L,   0 , 0, L^2];
    Ge(4,5:end) = [0     ,  3*L, -L^2, 0];
    Ge(5,6:end) = [-36   , -3*L,   0 ];
    Ge(6,7:end) = [0     , -3*L];
    Ge(7,8:end) = [-4*L^2];
    Ge = (Ge - Ge');
    Ge = Ge*jz/(30*L);  % Checked

    Ke(1,1:end) = [12   ,   0 , 0,  6*L, -12  ,   0 , 0, 6*L];
    Ke(2,2:end) = [12   , -6*L, 0,   0 , -12  , -6*L, 0];
    Ke(3,3:end) = [4*L^2,   0 , 0,  6*L, 2*L^2,   0];
    Ke(4,4:end) = [4*L^2, -6*L, 0,   0 , 2*L^2];
    Ke(5,5:end) = [12   ,   0 , 0, -6*L];
    Ke(6,6:end) = [12   ,  6*L, 0];
    Ke(7,7:end) = [4*L^2,   0 ];
    Ke(8,8:end) = [4*L^2];
    Ke = (Ke + Ke') - diag(diag(Ke));
    Ke = Ke*E*I/(L^3);  % Checked
   
    Fge = [0; m*g*L/2; -m*g*L^2/12; 0; 0; m*g*L/2; m*g*L^2/12; 0];
end