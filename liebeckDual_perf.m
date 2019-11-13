function analysis = liebeckDual_perf(inputs,design)

Wdot = inputs.Wdot;
V = inputs.V;
n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;
rho = inputs.rho;

xi = design.xi;
cF  = design.cF;
cR  = design.cR;
betaF = design.betaF;
betaR = design.betaR;
sigmaF = design.sigmaF;
sigmaR = design.sigmaR;

% [T, asound, patm, rho] = atmosisa(h);
R = D/2;
nu = 0.000014607;
rpm = 60*n;
omega = n*2*pi;
J = V/(n*D);
lambda = V/(omega*R);
S = pi*R^2; % Projected area of the helix

%% Polares
polares = 0; % 0 para não gerar polar, 1 para gerar
if (polares == 1)
    %     [pol foil] = xfoil('clark-y.dat',-5:0.1:15,mean(Re),0.3,'oper iter 300');
    %     for i = 1:length(x)
    pol = xfoil('clark-y.dat',-5:0.5:35,5e5,0.5,'oper iter 300');
    pol.CL_CD = pol.CL./pol.CD;
    %     end
    save('polares1.mat','pol');
else
    %     load('polares1.mat');
    load('polares3602.mat');
end


%% main loop

% V = 45:5:185;
V=153;
r = xi*R;

for i = 1:length(V)
    
    J(i) = V(i)/(n*D);
    omega = n*2*pi;
    lambda = V(i)/(omega*R);
    x = omega.*r/V(i);
    Pc1 = Wdot/(0.5*rho*V(i)^3*pi*R^2);
    Cp1 = Wdot/(rho*n^3*D^5);
        
    phi = atan(V(i)./(omega.*r));
    phiF = phi;
    phiR = phi;
    
    erro = inf;
    
    while erro > 1e-4
        phit = atan(tan(phi).*xi);
        x = omega*r/V(i);
        f = 0.5*B*(1-xi)./sin(phit);
        F = 2/pi*acos(exp(-f));
        G = F.*cos(phi).*sin(phi)*2;
        
        alphaF = betaF - phiF;
        alphaR = betaR - phiR;
        
        for j = 1:length(r)
            CLF(j) = interp1(pol.alpha,pol.CL,rad2deg(alphaF(j)));
            CDF(j) = interp1(pol.alpha,pol.CD,rad2deg(alphaF(j)));
            CLR(j) = interp1(pol.alpha,pol.CL,rad2deg(alphaR(j)));
            CDR(j) = interp1(pol.alpha,pol.CD,rad2deg(alphaR(j)));
        end
        
        eF = CDF./CLF;   
        eR = CDR./CLR;   
        CyF = CLF.*(cos(phiF) - eF.*sin(phiF));
        CxF = CLF.*(sin(phiF) + eF.*cos(phiF));
        CyR = CLR.*(cos(phiR) - eR.*sin(phiR));
        CxR = CLR.*(sin(phiR) + eR.*cos(phiR));
        
        KaF = CyF./(4*sin(phiF).^2);
        KbF = CxF./(4*cos(phiF).*sin(phiF));
        KaR = CyR./(4*sin(phiR).^2);
        KbR = CxR./(4*cos(phiR).*sin(phiR));
        
        aF = sigmaF.*KaF./(F - sigmaF.*KaF-sigmaR.*KaR);
        bF = sigmaF.*KbF./(F + 2*sigmaF.*KbF);
        aR = sigmaR.*KaR./(F - sigmaF.*KaF-sigmaR.*KaR);
%         bR=bF;
        bR = sigmaR.*KbR./(F + 2*sigmaR.*KbR);
%         aF = sigmaF.*KaF./(F - sigmaF.*KaF);
%         bF = sigmaF.*KbF./(F + sigmaF.*KbF);
%         aR = sigmaR.*KaR./(F - sigmaR.*KaR);
%         bR = sigmaR.*KbR./(F + sigmaR.*KbR);
%         aF(end) = 0;
%         bF(end) = 0;
%         aR(end) = 0;
%         bR(end) = 0;
        
        phinewF = atan(V(i).*(1+aF+aR)./(omega.*r.*(1-bF)));
        phinewR = atan(V(i).*(1+aF+aR)./(omega.*r.*(1+bF)));
        phinew = atan(V(i).*(1+aF+aR)./(omega.*r)); 

        erro = norm(abs(phiF-phinewF));
        
        phi = phinew;
        phiF = phinewF;
        phiR = phinewR;
        
    end
    
    WF = V(i).*(1+aF+aR)./sin(phiF);
    WR = V(i).*(1+aF+aR)./sin(phiR);
    
    dTF = 0.5*rho.*WF.^2*B/2.*cF.*CyF;
    TF(i) = trapz(r,dTF)
    dTR = 0.5*rho.*WR.^2*B/2.*cR.*CyR;
    TR(i) = trapz(r,dTR)
    T(i) = (TF(i) + TR(i))
    CtF(i) = TF(i)./(rho*n^2*D^4);
    CtR(i) = TR(i)./(rho*n^2*D^4);
    Ct(i) = T(i)./(rho*n^2*D^4);

    dQF = 0.5*rho.*WF.^2*B/2.*cF.*CxF.*r;
    QF = trapz(r,dQF);
    dQR = 0.5*rho.*WR.^2*B/2.*cR.*CxR.*r;
    QR = trapz(r,dQR);
    Q = (QF + QR);
    PF = omega*QF;
    PR = omega*QR;
    P = omega*Q;
    CpF(i) = PF./(rho*n^3*D^5);
    CpR(i) = PR./(rho*n^3*D^5);
    Cp(i) = P./(rho*n^3*D^5);
    etaF(i) = J(i)*CtF(i)/CpF(i);
    etaR(i) = J(i)*CtR(i)/CpR(i);
    eta(i) = J(i)*Ct(i)/Cp(i);

end


analysis.J = J;
analysis.CtF = CtF;
analysis.CtR = CtR;
analysis.Ct = Ct;
analysis.CpF = CpF;
analysis.CpR = CpR;
analysis.Cp = Cp;
analysis.etaF = etaF;
analysis.etaF = etaF;
analysis.eta = eta;
analysis.TF = TF;
analysis.TR = TR;
analysis.T = T;
