function [R, dRdU, dRda] = ROTBRNG_HBRESFUN(Ua, M, C, G, K, L, Fl, Fs, h, Nt, brngs, Cmat)
    Ua = [Cmat*Ua(1:end-1); Ua(end)];

    w = Ua(end);
    E = HARMONICSTIFFNESS(M, (C-w*G), K, w, h);
    dEdw = HARMONICSTIFFNESS(2*M/w, C/w-2*G, K*0, w, h);
    D1 = HARMONICSTIFFNESS(0, 1, 0, w, h);
    
    Nph = size(L,1);
    Nd = size(L,2);
    Nhc = 2*sum(h~=0) + sum(h==0);
    
    cst = TIMESERIES_DERIV(Nt, h, eye(Nhc), 0);
    sct = TIMESERIES_DERIV(Nt, h, D1, 1);
    
    props = struct('h0', brngs.h0(1), 'R', brngs.R(1), 'mu', brngs.mu(1), ...
        'L', brngs.L(1));
    
    Fnl = zeros(Nph*Nhc, 1);
    dFnldw = Fnl;
    Jnl = zeros(Nph*Nhc);
    for i=1:length(brngs.nds)
        ni = brngs.nds(i);
        exh = kron(eye(Nhc), L((ni-1)*4+1,:))*Ua(1:end-1);
        eyh = kron(eye(Nhc), L((ni-1)*4+2,:))*Ua(1:end-1);
        
        exyxdydt = [TIMESERIES_DERIV(Nt, h, [exh eyh], 0),...
            w*TIMESERIES_DERIV(Nt, h, [exh eyh], 1)];
       
        % Time Domain Force Calculation
        fxyt = zeros(Nt, 2);
        dfxydexyxdydt = zeros(2*Nt, 4);
        dfxydwt = fxyt;
        props = struct('h0', brngs.h0(i), 'R', brngs.R(i), 'mu', brngs.mu(i), ...
        'L', brngs.L(i));
        for n=1:Nt
            [fxyt(n,:), dfxydexyxdydt((n-1)*2+(1:2),:)] = ...
                BEARINGFORCES(exyxdydt(1,:), w, brngs.Nth, brngs.Nz, brngs.No, ... 
                props);
            dfxydwt(n,:) = dfxydexyxdydt((n-1)*2+(1:2),3:4)*(exyxdydt(n,3:4)')/w;
        end
        
        % Frequency Domain Transformation
        FXY = GETFOURIERCOEFF(h, fxyt);
        dFXdX = GETFOURIERCOEFF(h, dfxydexyxdydt(1:2:end, 1).*cst+dfxydexyxdydt(1:2:end, 3).*sct);
        dFYdX = GETFOURIERCOEFF(h, dfxydexyxdydt(2:2:end, 1).*cst+dfxydexyxdydt(2:2:end, 3).*sct);

        dFXdY = GETFOURIERCOEFF(h, dfxydexyxdydt(1:2:end, 2).*cst+dfxydexyxdydt(1:2:end, 4).*sct);
        dFYdY = GETFOURIERCOEFF(h, dfxydexyxdydt(2:2:end, 2).*cst+dfxydexyxdydt(2:2:end, 4).*sct);
        
        % Assembling into nonlinear forcing & jacobian
        Fnl(((ni-1)*4+1):Nph:end) = FXY(:,1);
        Fnl(((ni-1)*4+2):Nph:end) = FXY(:,2);
%         Fnl(((ni-1)*4+brngs.si(i))) = brngs.Fs(i);
        Jnl(((ni-1)*4+1):Nph:end, ((ni-1)*4+1):Nph:end) = dFXdX;
        Jnl(((ni-1)*4+1):Nph:end, ((ni-1)*4+2):Nph:end) = dFXdY;

        Jnl(((ni-1)*4+2):Nph:end, ((ni-1)*4+1):Nph:end) = dFYdX;
        Jnl(((ni-1)*4+2):Nph:end, ((ni-1)*4+2):Nph:end) = dFYdY;
        
        dFnldw(((ni-1)*4+1):Nph:end) = GETFOURIERCOEFF(h, dfxydwt(:,1));
        dFnldw(((ni-1)*4+2):Nph:end) = GETFOURIERCOEFF(h, dfxydwt(:,2));
    end
    
    R = E*Ua(1:end-1) + kron(eye(Nhc),L')*Fnl - w^2*Fl - Fs;
    dRdU = E + kron(eye(Nhc),L')*Jnl*kron(eye(Nhc),L);
    
    R = Cmat'*R;
    dRdU = Cmat'*dRdU*Cmat;
%     R = E*Ua(1:end-1) - kron(eye(Nhc),L')*Fnl - w^2*Fl;
%     dRdU = E - kron(eye(Nhc),L')*Jnl*kron(eye(Nhc),L);
    
    dRda = dEdw*Ua(1:end-1)+kron(eye(Nhc),L')*dFnldw-2*w*Fl;
    
    dRdX = [dRdU dRda];
end