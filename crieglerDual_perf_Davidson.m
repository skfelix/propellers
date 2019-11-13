function analysis = crieglerDual_perf_Davidson(inputs,design)

Wdot = inputs.Wdot;
V = inputs.V;
rho = inputs.rho;
n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;

x = design.xi;
bF  = design.cF;
bR  = design.cR;
betaF = design.betaF;
betaR = design.betaR;
sigmaF = design.sigmaF;
sigmaR = design.sigmaR;
CLF = design.CLF;
CLR = design.CLR;

crieglerPlots;

% [T, asound, patm, rho] = atmosisa(h);
R = D/2;
nu = 0.000014607;
rpm = 60*n;
omega = n*2*pi;
J = V/(n*D);
lambda = V/(omega*R);
S = pi*R^2; % Projected area of the helix
nP = 'dualRotation'; % 'singleRotation' or 'dualRotation'
interpMode = 'spline'; % 'linear' or 'spline

%% Polares
polares = 0; % 0 para não gerar polar, 1 para gerar
if (polares == 1)
    %     [pol foil] = xfoil('clark-y.dat',-5:0.1:15,mean(Re),0.3,'oper iter 300');
    %     for i = 1:length(x)
    pol = xfoil('clark-y.dat',-5:0.5:35,5e5,0.3,'oper iter 300')
    pol.CL_CD = pol.CL./pol.CD;
    %     end
    save('polares1.mat','pol');
else
    %     load('polares1.mat');
    load('polares3602.mat');
end

%%

V = (50:5:185)*.305;
% V = 153;
for jj = 1:length(V)
    J(jj) = V(jj)/(n*D);
    
    W = sqrt(V(jj)^2 + (omega.*x.*R).^2);
    nu = 0.000014607;
    Re = W.*bF/nu;
    Mach = W/340;
    
    %% Induced Velocity
    Pct = Wdot/(0.5*rho*V(jj)^3*pi*(D/2)^2);
    for wbar = 0:0.0001:0.5
        Jw(jj) = J(jj)*(1+wbar);
        [k,ek] = getKappa(nP, interpMode, B, Jw(jj));
        Pc = 2*k*wbar*(1+wbar)*(1+ek*wbar);
        if abs(Pct - Pc) < 1e-3
            wout(jj) = wbar;
            break;
        end
    end
    
    phi = atan(1/pi*J(jj)*(1+0.5*wbar)./x);
    phi0 = atan(J(jj)./(pi.*x));
    x1 = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
    [KxVec1] = getKx(nP, interpMode, B, Jw(jj)); % KxVec comes for x1 = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
    KxVec = interp1(x1,KxVec1,x,'linear', 'extrap');
    
    %% Tip Loss Correction
    
    p = 0.5.*sin(phi0).*cos(phi0);
    q = 0.5./sin(phi0);
    r = cos(phi0).^2 - sin(phi0).^2;
    s = Jw(jj)./(pi*x).*(KxVec.*sin(phi0));
    
    
    X0 = q.*s./(2.*p - q.*r.*s);
    b_line = 1./(4*X0.*sin(phi0));
    
    erro = inf;
    %%
    while erro > 1e-4
        
        alphaiF = b_line.*sigmaF.*CLF.*(1+X0.*cos(phi0).^2);
        alphaiR = b_line.*sigmaR.*CLR.*(1+X0.*(cos(phi0).^2 - 2*sin(phi0).^2));
        alphaF = betaF - phi0 - alphaiF;
        alphaR = betaR - phi0 - alphaiR;
        phiF = betaF - alphaF;
        phiR = betaR - alphaR;
        
        for i = 1:length(x)
            CLFnew(i) = interp1(pol.alpha,pol.CL,rad2deg(alphaF(i)));
            CLRnew(i) = interp1(pol.alpha,pol.CL,rad2deg(alphaR(i)));
            CDF(i) = interp1(pol.alpha,pol.CD,rad2deg(alphaF(i)));
            CDR(i) = interp1(pol.alpha,pol.CD,rad2deg(alphaR(i)));
        end
        
        erro = abs(CLFnew(i) - CLF(i));
        
        CLF = CLFnew;
        CLR = CLRnew;
    end
    %%
    eF = CDF./CLF;
    eR = CDR./CLR;
    phiq0F = sin(phi0) + eF.*cos(phi0);
    phiq0R = sin(phi0) + eR.*cos(phi0);
    
    sigmaCLF = sigmaF.*CLF;
    sigmaCLR = sigmaR.*CLR;
    
    dCpdxF = pi^4/4*sec(phi0).^2.*x.^4.*sigmaCLF.*phiq0F;
    dCpdxR = pi^4/4*sec(phi0).^2.*x.^4.*sigmaCLR.*phiq0R;
%     CpF(jj) = trapz(x,dCpdxF);
%     CpR(jj) = trapz(x,dCpdxR);
%     Cp(jj) = CpF(jj) + CpR(jj);

    
%     dQ = 0.5*rho.*W.^2*B.*bF.*CLF.*(sin(phi) + eF.*cos(phi0)).*(x*R);
%     dQ = 0.5*rho.*W.^2*B.*bF.*0.5.*(CLF+CLR).*(sin(phi) + eF.*cos(phi)).*(x*R);
    dQ = pi*rho*(x*R).^2.*sigmaCLF.*W.^2.*(sin(phiF) + eF.*cos(phiF)) + ...
         pi*rho*(x*R).^2.*sigmaCLR.*W.^2.*(sin(phiR) + eR.*cos(phiR));
    Q = trapz(x*R,dQ);
    P = omega*Q;
    Cp(jj) = P./(rho*n^3*D^5);
    
%         dT = 0.5*rho.*W.^2*B.*bF.*CLF.*(cos(phi0) - eF.*sin(phi0));
    %     dT = 0.5*rho.*W.^2*B.*bF.*0.5.*(CLF+CLR).*(cos(phi) - eF.*sin(phi));
    dT =  pi*rho*(x*R).*sigmaCLF.*W.^2.*(cos(phiF) - eF.*sin(phiF))+ ...
          pi*rho*(x*R).*sigmaCLR.*W.^2.*(cos(phiR) - eR.*sin(phiR));
    
    T(jj) = trapz(x*R,dT);
    Ct(jj) = T(jj)./(rho*n^2*D^4);
    
    eta(jj) = J(jj)*Ct(jj)/Cp(jj);
    
    %     CtF(jj) = etaF(jj)*CpF(jj)/J(jj);
    %     CtR(jj) = etaR(jj)*CpR(jj)/J(jj);
    %     Ct(jj) = eta(jj)*Cp(jj)/J(jj);
    %     TR(jj) = CtR(jj)*(rho*n^2*D^4);
    %     TF(jj) = CtF(jj)*(rho*n^2*D^4);
    %     T(jj) = Ct(jj)*(rho*n^2*D^4);
    %
end


analysis.J = J;
% analysis.CtF = CtF;
% analysis.CtR = CtR;
analysis.Ct = Ct;
analysis.CpF = Cp;
% analysis.CpR = CpR;
analysis.Cp = Cp;
analysis.eta = eta;
% analysis.etaF = etaF;
% analysis.etaR = etaR;
% analysis.TF = TF;
% analysis.TR = TR;
analysis.T = T;
% analysis.a1 = a1;