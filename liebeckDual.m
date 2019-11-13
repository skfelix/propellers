% function out = liebeckDual(inputs,polares,propType,CL_design,xi)

Wdot = inputs.Wdot;
V = inputs.V;
n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;
rho = inputs.rho;

% [T, asound, patm, rho] = atmosisa(h);
R = D/2;
nu = 0.000014607;
rpm = 60*n;
omega = n*2*pi;
J = V/(n*D);
lambda = V/(omega*R);
S = pi*R^2; % Projected area of the helix

Cp = Wdot/(rho*n^3*D^5);
Cq = Cp/omega;
Pc = Wdot/(0.5*rho*V^3*pi*R^2);
%% Polares - Executar somente se mudar as condições de contorno
% polares = 0; % 0 para não gerar polar, 1 para gerar
if (polares == 1)
    [p f] = xfoil('clark-y.dat',0:0.05:8,5e5,0.25,'oper iter 100');
    p.CL_CD = p.CL./p.CD;
    save('polar_clarky.mat','p');
else
    load('polar_clarky.mat');
end

%% Main Loop
r = xi*R;

zeta = 0.2;
erro = inf;

while erro > 1e-4
    %% Prandtl Factor
    phit = atan(lambda*(1+0.5*zeta));
    
    x = omega*r/V;
    f = 0.5*B*(1-xi)/sin(phit);
    F = 2/pi*acos(exp(-f));
    phi = atan(tan(phit)./xi);
%     Fh = 2/pi*acos(exp(-B/2*(x-0.2)./(0.2*sin(phi))));
%     F = F.*Fh;
    g= 0.83;
%     g = F.*sin(phi);%.*cphi)*2;
% g=[0.8212    0.8195    0.8174    0.8149    0.8116    0.8077    0.8027    0.7965    0.7889    0.7794    0.7675    0.7528    0.7344    0.7114    0.6825    0.6460    0.5993    0.5385    0.4562    0.3346         0];
%      g = F.*sin(phi).*cos(phi)*2;
%      g = .82./F;
     
    
    %% Prop Type
    
    if propType == 1 % max CL/CD 
        for i = 1:length(r)
            cl_cd_max = max(p.CL_CD);
            alpha(i) = interp1(p.CL_CD,p.alpha,cl_cd_max);
            CL(i) = interp1(p.alpha,p.CL,alpha(i));
            CD(i) = interp1(p.alpha,p.CD,alpha(i));
        end
        e = 1/cl_cd_max;
    elseif propType == 2 % min CD
        for i = 1:length(r)
            cd_min = min(p.CD);
            alpha(i) = interp1(p.CD,p.alpha,cd_min);
            CL(i) = interp1(p.alpha,p.CL,alpha(i));
            CD(i) = cd_min;
        end
        e = CD./CL;
    elseif propType == 3 % fixed CL | CL design
%         CL_design = 0.7;
        for i = 1:length(r)
            alpha(i) = interp1(p.CL,p.alpha,CL_design);
            CL(i) = CL_design;
            CD(i) = interp1(p.alpha,p.CD,alpha(i));
        end
        e = CD./CL;
    elseif propType == 4 % Takeoff
        CL_takeoff = 0.5;
        for i = 1:length(r)
            alpha(i) = interp1(p.CL,p.alpha,CL_takeoff);
            CL(i) = CL_takeoff;
            CD(i) = interp1(p.alpha,p.CD,alpha(i));
        end
        e = CD./CL;
    end
        
    Cy = CL.*(cos(phi) - e.*sin(phi));
    Cx = CL.*(sin(phi) + e.*cos(phi));
    %%
    
    WcF = 4*pi.*g.*F.*V.^2*zeta./(omega*CL*B/2);
    WcR = 4*pi.*g.*F.*V.^2*zeta./(omega*CL*B/2);
    ReF = WcF/nu;
    ReR = WcR/nu;
    Re = (WcF+WcR)/nu;
    
    bF = 0.5*zeta.*g./x.^2; 
    bR = bF;
    z = 0.5*zeta.*g./x.^2; 
    aR = (-2 + sqrt(4+8.*(1-z).*zeta.*g.*(1+(1-z)./(1+z))))./ ...
            (4*(1+(1-z)./(1+z)));
    aF = aR.*(1-z)./(1+z);
%     aF(end)=0;aR(end)=0;
%     bF(end)=0;bR(end)=0;
    
    phiF = atan(V*(1+aF+aR)./(omega*r.*(1-bF)));
    phiR = atan(V*(1+aF+aR)./(omega*r.*(1+bR)));
    betaF = deg2rad(alpha)+phiF;
    betaR = deg2rad(alpha)+phiR;
    phiw = atan(V*(1+aF+aR)./(omega*r));
    
    WF = V*(1+aF+aR)./sin(phiF);
    WR = V*(1+aF+aR)./sin(phiR);
    cF = WcF./WF;
    cR = WcR./WR;
    
    dI1 = 4*xi.*g.*F;
    dI2 = 2*xi.*g.^2.*F./(x.^2); %% Conferir
    dJ1 = 4*xi.*g.*F;
    dJ2 = 4*xi.*g.^2.*F./(x.*tan(phiw)); %% Conferir .*tan
    
    I1 = trapz(xi,dI1);
    I2 = trapz(xi,dI2);
    J1 = trapz(xi,dJ1);
    J2 = trapz(xi,dJ2);
    zetanew = (0.5*J1./J2)*(sqrt(1+4*(Pc/2)*J2./J1.^2)-1);
    
    TcF = I1.*zeta - I2.*zeta.^2;
    TcR = I1.*zeta + I2.*zeta.^2;
    Tc = TcF+TcR;
    
    erro = abs(zeta - zetanew);
    
    zeta = zetanew;
    
end

cCLF = cF.*CL;
cCLR = cR.*CL;
sigmaF = (B/2)*cF./(2*pi.*r);
sigmaR = (B/2)*cR./(2*pi.*r);
sigmaCLF = sigmaF.*CL;
sigmaCLR = sigmaR.*CL;
TF = TcF*0.5*rho*V^2*S;
TR = TcR*0.5*rho*V^2*S;
Ts = Tc*0.5*rho*V^2*S;
Ct = Ts./(rho*n^2*D^4);
CtF = TF./(rho*n^2*D^4);
CtR = TR./(rho*n^2*D^4);
CpF = Cp/2;
CpR = Cp/2;
etaF = J*CtF/CpF;
etaR = J*CtR/CpR;
eta = J*Ct/Cp;

dQF = 0.5*rho.*WF.^2*B/2.*cF.*Cx.*r;
QF = trapz(r,dQF);
dQR = 0.5*rho.*WR.^2*B/2.*cR.*Cx.*r;
QR = trapz(r,dQR);
Q = QF + QR;
PF = omega*QF;
PR = omega*QR;
P = omega*Q;
Cp2 = P./(rho*n^3*D^5);

dTF = 0.5*rho.*WF.^2*B/2.*cF.*Cy;
TF = trapz(r,dTF);
dTR = 0.5*rho.*WR.^2*B/2.*cR.*Cy;
TR = trapz(r,dTR);
T = TF + TR;
Ct2 = T./(rho*n^2*D^4);

eta = J*Ct2/Cp2;



beta75F = interp1(xi,betaF,0.75);
beta75  = interp1(xi,deg2rad(alpha)+phiw,0.75);
betadeg75 = rad2deg(beta75);
beta75R = interp1(xi,betaR,0.75);
pitch = 2*pi*(xi*R).*tan(deg2rad(alpha)+phiw);
pitch75 = 2*pi*(0.75*R)*tan(beta75);

solidityF = B/2*interp1(xi,cF,0.75)/(2*pi*0.75*R);
solidityR = B/2*interp1(xi,cR,0.75)/(2*pi*0.75*R);
solidity = solidityF + solidityR;
TAFF = B/2*1e5/(32*R^5)*trapz(xi*R, cF.*(xi*R).^3);
TAFR = B/2*1e5/(32*R^5)*trapz(xi*R, cR.*(xi*R).^3);
TAF = TAFF + TAFR;

%%
% deltapa = trapz(xi*D,0.5*rho*((V*(1+aF+aR)).^2-V^2))
% deltapr = trapz(r,rho*(omega-0.5*bF*omega).*bF.*omega.*r.^2 + rho*(omega+0.5*bR*(omega)).*bR.*(omega).*r.^2)
% % deltap = trapz(r,0.5*rho*omega^2.*r.^2.*(bF+bR -0.5*(bF.^2 -bR.^2)))
% trapz(r,rho.*omega^2.*r.^2.*(bF+bR))
% deltap = deltapa - deltapr
% % TaF = trapz(r,4*rho*pi*r.*V^2.*aF.*(1+aF+aR).*F)
% % TaR = trapz(r,4*rho*pi*r.*V^2.*aR.*(1+aF+aR).*F)
% % Ta = TaF + TaR
% % 
% % % Ta = trapz(r,4*rho*pi*r.*V^2.*a.*(1+a).*F)
% trapz(r,rho*pi*r.*V^2.*aF.*(2+3*aF)) + deltap*pi*R^2 + ...
%     trapz(r,rho*pi*r.*(V*(1+aF)).^2.*aR.*(2+3*aR))
% % trapz(r,rho*pi*r.*V^2.*((1+2*aF+2*aR).^2 -(1+aF+aR).^2).*F) + deltap*pi*R^2 %+ ...
% %     trapz(r,rho*pi*r.*V^2.*((1+2*aF+2*aR).^2 -(1+aF+aR).^2))
%%

out.pitch = pitch;
out.pitch75 = pitch75;
out.xi = xi;
out.cF = cF;
out.cR = cR;
out.betaF = betaF;
out.betaR = betaR;
out.betadegF = rad2deg(betaF);
out.betadegR = rad2deg(betaR);
out.beta75F = beta75F;
out.beta75R = beta75R;
out.beta75 = beta75;
out.betadeg75 = betadeg75;
out.sigmaF = sigmaF;
out.sigmaR = sigmaR;
out.sigmaCLF = sigmaCLF;
out.sigmaCLR = sigmaCLR;
out.cCLF = cCLF;
out.cCLR = cCLR;
out.CLF = CL;
out.CLR = CL;
out.CtF = CtF;
out.CtR = CtR;
out.Ct = Ct;
out.CpF = CpF;
out.CpR = CpR;
out.Cp = Cp;
out.CqF = CpF/omega;
out.CqR = CpR/omega;
out.Cq = Cp/omega;
out.TF = TF;
out.TR = TR;
out.T = T;
% out.PF = PF;
% out.PR = PR;
% out.P = P;
% out.QF = QF;
% out.QR = QR;
% out.Q = Q;
out.J = J;
out.etaF = etaF;
out.etaR = etaR;
out.eta = eta;
out.TAFF = TAFF;
out.TAFR = TAFR;
out.TAF = TAF;
out.solidityF = solidityF;
out.solidityR = solidityR;
out.solidity = solidity;
