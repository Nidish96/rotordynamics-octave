function [Me, Ke] = EBBEAM3D_MATS(rho, E, A, I1, I2, L)
%EBBEAM3D_MATS returns the element matrices
%
% USAGE:
%   [Me, Ke] = EBBEAM3D_MATS(rho, E, A, I1, I2, L);
% INPUTS:
%   rho     : density
%   E       : Young's modulus
%   A       : Cross sectional area
%   I1,I2   : Second moment of cross section (about bending axes)
%   L       : Element length
% OUTPUTS:
%   Me      : 10x10 Element Mass matrix
%   Ke      : 10x10 Element Stiffness matrix

    Me = sparse(10,10);
    Ke = sparse(10,10);

    % Mass Matrix
    Me([1 6], [1 6]) = rho*A*L/6*[2 1; 1 2];  % Axial
    Me([2 3 7 8], [2 3 7 8]) = rho*A*L/420*...
        [156 22*L 54 -13*L;
        22*L 4*L^2 13*L -3*L^2;
        54 13*L 156 -22*L;
        -13*L -3*L^2 -22*L 4*L^2] + rho*I1/(30*L)*...
        [36 3*L -36 3*L;
        3*L 4*L^2 -3*L -L^2;
        -36 -3*L 36 -3*L;
        3*L -L^2 -3*L 4*L^2];  % Bending 1
    Me([4 5 9 10], [4 5 9 10]) = rho*A*L/420*...
        [156 22*L 54 -13*L;
        22*L 4*L^2 13*L -3*L^2;
        54 13*L 156 -22*L;
        -13*L -3*L^2 -22*L 4*L^2] + rho*I2/(30*L)*...
        [36 3*L -36 3*L;
        3*L 4*L^2 -3*L -L^2;
        -36 -3*L 36 -3*L;
        3*L -L^2 -3*L 4*L^2];  % Bending 2
    
    % Stiffness Matrix
    Ke([1 6], [1 6]) = A*E/L*[1 -1;-1 1];  % Axial
    Ke([2 3 7 8], [2 3 7 8]) = 2*E*I1/L^3*...
        [6   3*L   -6   3*L;
         3*L 2*L^2 -3*L L^2;
         -6  -3*L  6    -3*L;
         3*L L^2   -3*L 2*L^2];  % Bending 1
    Ke([4 5 9 10], [4 5 9 10]) = 2*E*I2/L^3*...
        [6 3*L -6 3*L;
        3*L 2*L^2 -3*L L^2;
        -6 -3*L 6 -3*L;
        3*L L^2 -3*L 2*L^2];  % Bending 2

end