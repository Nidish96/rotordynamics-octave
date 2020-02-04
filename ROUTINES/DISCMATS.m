function [M, G, FC, FS, MT, MR] = DISCMATS(m, J, Jz, Jxz, Jyz, ax, ay)
    MT = diag([1, 1, 0, 0])*m;
    MR = diag([0, 0, 1, 1])*J;
    M = MT+MR;
    G = zeros(4,4);
    G(3:4, 3:4) = [0 -1;1 0]*Jz;
    FC = [m*ax; m*ay; Jxz; Jyz];
    FS = [-m*ay; m*ax; -Jyz; Jxz];
end