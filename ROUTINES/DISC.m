classdef DISC
    properties
        nd  % Node
        Rs  % Shaft Radius
        Rd  % Disc Radius
        t  % Thickness
        rho  % Density
        ax = 0  % X offset
        ay = 0  % Y offset 
        m  % Mass
        g
        Jxx
        Jyy
        Jzz
        Jxz = 0
        Jyz = 0
    end
    methods
        function obj = DISC(nd, Rs, Rd, t, rho, ax, ay, g)
            obj.nd = nd;
            obj.Rs = Rs;
            obj.Rd = Rd;
            obj.t  = t;
            obj.rho= rho;
            obj.ax = ax;
            obj.ay = ay;
            obj.g  = g;
            
            obj.m   = pi*rho*t*(Rd^2-Rs^2);
            obj.Jzz = obj.m*(Rd^2+Rs^2)/2 + obj.m*(ax^2+ay^2);  % Polar M.I.
            obj.Jxx = obj.m*(Rd^2+Rs^2)/4 + obj.m*ay*ay;  % Dimetral x
            obj.Jyy = obj.m*(Rd^2+Rs^2)/4 + obj.m*ax*ax;  % Dimetral y
            obj.Jxz = 0;
            obj.Jyz = 0;
        end
        function [M, G, FC, FS, Fg, MT, MR] = MATS(obj)
            MT = diag([1, 1, 0, 0])*obj.m;
            MR = diag([0, 0, obj.Jxx, obj.Jyy]);
            M = MT+MR;
            G = zeros(4,4);
            G(3:4, 3:4) = [0 -1;1 0]*obj.Jzz;
            FC = [obj.m*obj.ax; obj.m*obj.ay; obj.Jxz; obj.Jyz];
            FS = [-obj.m*obj.ay; obj.m*obj.ax; -obj.Jyz; obj.Jxz];
            Fg = [0; obj.m*obj.g; 0; 0];
        end
    end
end
