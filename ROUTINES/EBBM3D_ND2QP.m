function [Q, T] = EBBM3D_ND2QP(Les, No)
%EBBM3D_ND2QP returns nodal to quadrature interpolation matrix and qp to
%nodal integration matrix
%
% USAGE:
% INPUTS:
% OUTPUTS:

    [xi, w] = LGWT(No, -1, 1);
    [xi, si] = sort(xi);
    w = w(si);
    Ne = length(Les);
    Q  = zeros(Ne*No*3, (Ne+1)*5);
    T  = zeros((Ne+1)*5, Ne*No*3);
    for e=1:Ne
        NS = SFUNCS(xi, Les(e));
        Q((e-1)*(No*3)+(1:No*3), (e-1)*5+(1:10)) = NS;
        
        T((e-1)*5+(1:10), (e-1)*No*3+(1:No*3)) = T((e-1)*5+(1:10), (e-1)*No*3+(1:No*3)) + ...
            (kron(w, ones(3,1)).*NS*Les(e)/2)';
    end
end

function [Ns] = SFUNCS(xi, Le)
    Np = length(xi);
    Ns = zeros(Np*3, 10);
    Ns(1:3:end, [1 6])     = [1-xi, 1+xi]/2;
    Ns(2:3:end, [2 3 7 8]) = [2*(xi-1).^2.*(xi+2), Le*(xi-1).^2.*(xi+1),...
        -2*(xi-2).*(xi+1).^2, Le*(xi-1).*(xi+1).^2]/8;
    Ns(3:3:end, [4 5 9 10]) = [2*(xi-1).^2.*(xi+2), Le*(xi-1).^2.*(xi+1),...
        -2*(xi-2).*(xi+1).^2, Le*(xi-1).*(xi+1).^2]/8;
end