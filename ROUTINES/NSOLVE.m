function [U, r, eflag, it, jc] = NSOLVE(func, u0, opt)
%NSOLVE Uses Newton iterations to solve
%
% USAGE:
% INPUTS:
% OUTPUTS:

    [r0, j0] = func(u0);
    du0 = -j0\r0;
    e0  = abs(r0'*du0);
    if (e0 < eps)
        e0 = 1.0;
        r0 = 1.0;
    end
    r   = r0;
    e   = e0;
    u   = u0;
    du  = du0;
    it  = 0;
    while ((e/e0 > opt.tole && vecnorm(r)/vecnorm(r0) > opt.tolr) && it<=opt.ITMAX)
        u  = u + du;
        it = it+1;
        
        [r, jc] = func(u);
        du = -jc\r;
        e  = abs(r'*du);
        fprintf('%d %e %e\n', it, e/e0, vecnorm(r)/vecnorm(r0));
    end
	U = u;
    R = r;
    if (e/e0 > opt.tole)
        eflag = -1;
    else
        eflag = 1;
    end
end