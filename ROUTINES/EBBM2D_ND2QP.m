function [Q, T] = EBBM2D_ND2QP(Les, No)
%EBBM2D_ND2QP returns nodal to quadrature interpolation matrix and qp to
%nodal integration matrix
%
% USAGE:
% INPUTS:
% OUTPUTS:

    [xi, w] = LGWT(No, -1, 1);
    [xi, si] = sort(xi);
    w = w(si);
    Ne = length(Les);
    Q  = zeros(Ne*No*2, (Ne+1)*3);
    T  = zeros((Ne+1)*3, Ne*No*2);
    for e=1:Ne
        NS = SFUNCS(xi, Les(e));
        Q((e-1)*(No*2)+(1:No*2), (e-1)*3+(1:6)) = NS;
        
        T((e-1)*3+(1:6), (e-1)*No*2+(1:No*2)) = T((e-1)*3+(1:6), (e-1)*No*2+(1:No*2)) + ...
            (kron(w, ones(2,1)).*NS*Les(e)/2)';
    end
end

function [Ns] = SFUNCS(xi, Le)
    Np = length(xi);
    Ns = zeros(Np*2, 6);
    Ns(1:2:end, [1 4])     = [1-xi, 1+xi]/2;
    Ns(2:2:end, [2 3 5 6]) = [2*(xi-1).^2.*(xi+2), Le*(xi-1).^2.*(xi+1),...
                            -2*(xi-2).*(xi+1).^2, Le*(xi-1).*(xi+1).^2]/8;
end