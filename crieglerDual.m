function out = crieglerDual(inputs,polares,propType,CL_design,xi)

Wdot = inputs.Wdot;
V = inputs.V;
rho = inputs.rho;
n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;

crieglerPlots;

% [T, asound, patm, rho] = atmosisa(h);
nu = 0.000014607;
rpm = 60*n;
J = V/(n*D);
nP = 'dualRotation'; % 'singleRotation' or 'dualRotation'
interpMode = 'spline'; % 'linear' or 'spline
R = D/2;
J = V/(n*D);
omega = n*2*pi;
F = pi*R^2; % Projected area of the helix
W = sqrt(V^2 + (omega*R)^2);

Cpt = Wdot/(rho*n^3*D^5);
Pct = Wdot/(0.5*rho*V^3*pi*(D/2)^2);
%% Interference velocity

for wbar = 0:0.0001:1
    Jw = J*(1+wbar);
    [k,ek] = getKappa(nP, interpMode, B, Jw);
    Pc = 2*k*wbar*(1+wbar)*(1+ek*wbar);
    if abs(0.99*Pct - Pc) < 1e-4
        break;
    end
end

Pc_k = Pc/k;
eta = getEff(Pc_k, ek, interpMode);

%% Tables
x = xi;
w = wbar*V;
x_Kx = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
[KxVec1] = getKx(nP, interpMode, B, Jw); 
KxVec = interp1(x_Kx,KxVec1,x,'linear', 'extrap');
phi =  atan(J./(pi.*x).*(1+0.5*wbar));
phi0 = atan(J./(pi.*x));
phiF = atan(J./(pi.*x).*(1+0.5*wbar*(1+0.5*k.*(tan(phi).^2))));
phiR = atan(J./(pi.*x).*(1+0.5*wbar*(1-0.5*k.*(tan(phi).^2))));
sigmaCLF = J./(pi.*x).*(1+wbar).*wbar.*sin(phi0)./(1+0.25*k*wbar.*sin(phi0).^2).*KxVec;
sigmaCLR = J./(pi.*x).*(1+wbar).*wbar.*sin(phi0)./(1+0.75*k*wbar.*sin(phi0).^2).*KxVec;
bCLF = sigmaCLF.*2.*pi.*(x.*R)./(B/2);
bCLR = sigmaCLR.*2.*pi.*(x.*R)./(B/2);

nu = 0.000014607;
vr = sqrt(V^2 + (omega.*x.*R).^2);
% Re = vr.*b_guess/nu;
Mach = vr/340;

%% Polares - Executar somente se mudar as condições de contorno
% polares = 0; % 0 para não gerar par, 1 para gerar
if (polares == 1)
    for i = 1:length(x)
        [p f] = xfoil('clark-y.dat',0:0.2:10,5e5,0.3,'oper iter 300');
        p.CL_CD = p.CL./p.CD;
        p = p;
        foilRe{i} = f;
    end
%     save('polar_clarky.mat','pRe','foilRe');
else
    load('polar_clarky.mat');
end

%% Prop Type

if propType == 1 % max CL/CD
    for i = 1:length(x)
        idx = find(p.CL_CD==max(p.CL_CD));
        alpha(i) = p.alpha(idx); idx=idx(1);
        CLF(i) = p.CL(idx);
        CDF(i) = p.CD(idx);
        CLR(i) = p.CL(idx);
        CDR(i) = p.CD(idx);
    end
    e = 1/p.CL_CD(idx);
    bF = bCLF./CLF;
    bR = bCLR./CLR;
elseif propType == 2 % min CD
    for i = 1:length(x)
        idx = find(p.CD==min(p.CD)); idx=idx(1);
        alpha(i) = p.alpha(idx);
        CLF(i) = p.CL(idx);
        CDF(i) = p.CD(idx);
        CLR(i) = p.CL(idx);
        CDR(i) = p.CD(idx);
    end
    e = 1/p.CL_CD(idx);
    bF = bCLF./CLF;
    bR = bCLR./CLR;
elseif propType == 3 % fixed CL | CL design
    for i = 1:length(x)
        alpha(i) = interp1(p.CL,p.alpha,CL_design);
        CLF(i) = CL_design;
        CDF(i) = interp1(p.alpha,p.CD,alpha(i));
        CLR = CLF;
        CDR = CDF;
    end
    bF = bCLF./CLF;
    bR = bCLR./CLR;
elseif propType == 4 % Takeoff
    for i = 1:length(x)
        alpha(i) = interp1(p.CL,p.alpha,CL_design);
        CLF(i) = CL_design;
        CDF(i) = interp1(p.alpha,p.CD,alpha(i));
        CLR = CLF;
        CDR = CDF;
    end
    bF = bCLF./CLF;
    bR = bCLR./CLR;
end

%%
betaF = phiF + deg2rad(alpha);
betaR = phiR + deg2rad(alpha);
sigmaF = sigmaCLF./CLF;
sigmaR = sigmaCLR./CLR;

lambda_g = J/pi;
sigmaCDF = sigmaCLF.*CDF./CLF;
sigmaCDR = sigmaCLR.*CDR./CLR;
trF = 2/lambda_g^2*trapz(x,sigmaCDF.*x.^3./sin(phiF));
trR = 2/lambda_g^2*trapz(x,sigmaCDR.*x.^3./sin(phiR));
taF = 2*trapz(x,sigmaCDF.*x./sin(phiF));
taR = 2*trapz(x,sigmaCDR.*x./sin(phiR));
cs = 2*k*wbar*(1 + wbar*(0.5 + ek)); % total thrust coeff
csT = cs - taF - taR; % net thrust coeff
Pci = Pc + trF + trR; % total power coeff
eta_D = csT/Pci;
% TF = rho*D/4*V^3*wbar*(1+wbar)/n.* ...
%     trapz(x,(1+1/4.*k.*wbar.*sin(phi0).^2).*cos(phiF)./sin(phi0).*KxVec); % total thrust
% TR = rho*D/4*V^3*wbar*(1+wbar)/n.* ...
%     trapz(x,(1+3/4.*k.*wbar.*sin(phi0).^2).*cos(phiR)./sin(phi0).*KxVec); % total thrust
% T = TF+TR;
Ti = 1/2*rho*V^2*F*eta*Pci;
Ts = 1/2*rho*V^2*F*csT;
Ct = Ts/(rho*n^2*D^4);
P = 1/2*rho*V^3*F*Pci;
Cp = P/(rho*n^3*D^5);
Cq = Cp/(2*pi);
Q = Cq*rho*n^2*D^5;
solidityF = B/2*interp1(x,bF,0.75)/(2*pi*0.75*R);
solidityR = B/2*interp1(x,bR,0.75)/(2*pi*0.75*R);
solidity = solidityF+solidityR;
TAFF = B/2*1e5/(32*R^5)*trapz(x*R, bF.*(x*R).^3);
TAFR = B/2*1e5/(32*R^5)*trapz(x*R, bR.*(x*R).^3);
TAF = TAFF + TAFR;
beta75F = interp1(x,betaF,0.75);
beta75  = interp1(x,deg2rad(alpha)+phi,0.75);
betadeg75 = rad2deg(beta75);
beta75R = interp1(x,betaR,0.75);

pitch = 2*pi*(xi*R).*tan(deg2rad(alpha)+phi);
pitch75 = 2*pi*(0.75*R)*tan(beta75);

% e = atan(CDF./CLF);
% dCqF_dx = pi/8*J^2*(1+1/4*k*wbar*sin(phi0).^2).^2.*sigmaCLF.*cos(phiF)./(sin(phi0).^2).*(tan(phiF) + tan(e)).*x.^2;
% CqF = trapz(x,dCqF_dx);
% dCqR_dx = pi/8*J^2*(1+3/4*k*wbar*sin(phi0).^2).^2.*sigmaCLR.*cos(phiR)./(sin(phi0).^2).*(tan(phiR) + tan(e)).*x.^2;
% CqR = trapz(x,dCqR_dx);
% % dCtF_dx = pi/4*J^2*((1+1/4*k*wbar*sin(phi0).^2)./sin(phi0)).^2.*sigmaCLF.*(cos(phi) - e.*sin(phi)).*x;
% dCtF_dx = 2./x.*(1 - tan(phiF).*tan(e))./(tan(phiF) + tan(e)).*dCqF_dx;
% CtF = trapz(x,dCtF_dx);
% % dCtR_dx = pi/4*J^2*((1+3/4*k*wbar*sin(phi0).^2)./sin(phi0)).^2.*sigmaCLF.*(cos(phi) - e.*sin(phi)).*x;
% dCtR_dx = 2./x.*(1 - tan(phiR).*tan(e))./(tan(phiR) + tan(e)).*dCqR_dx;
% CtR = trapz(x,dCtR_dx);
% CpF = CqF*2*pi;
% CpR = CqR*2*pi;
% Cq = CqF + CqR;
% Ct = CtF + CtR;
% Cp = Cq*2*pi;
% eta = J*Ct/Cp;
% etaF = J*CtF/CpF;
% etaR = J*CtR/CpR;
% T = Ct*(rho*n^2*D^4);
% QF = CqF*(rho*n^2*D^5);
% QR = CqR*(rho*n^2*D^5);
% Q = QF + QR;
% PF = CpF*(rho*n^3*D^5);
% PR = CpR*(rho*n^3*D^5);
% P = PF + PR;

% eF = atan(CDF./CLF);
% eR = atan(CDR./CLR);
% dQF = pi*rho*(x*R).^2.*sigmaCLF.*W.^2.*(sin(phiF) + eF.*cos(phiF));
% dQR = pi*rho*(x*R).^2.*sigmaCLR.*W.^2.*(sin(phiR) + eR.*cos(phiR));
% QF = trapz((x*R),dQF);
% QR = trapz((x*R),dQR);
% Q = QF + QR;
% PF = omega*QF;
% PR = omega*QR;
% P = PR+PF;
% Cp = P./(rho*n^3*D^5);
% 
% dTF =  pi*rho*(x*R).*sigmaCLF.*W.^2.*(cos(phiF) - eF.*sin(phiF));
% dTR =  pi*rho*(x*R).*sigmaCLR.*W.^2.*(cos(phiR) - eR.*sin(phiR));
% TF = trapz((x*R),dTF);
% TR = trapz((x*R),dTR);
% T = TF + TR;
% 
% Ct = T./(rho*n^2*D^4);
% 
% eta = J*Ct/Cp;

% fprintf('=================== Propeller Characteristics ===================\n')
% fprintf('x\t\tK(x)\ttan(phiF)\ttan(phiR)\tsigmaCLF\tsigmaCLR\tbCLF\tbCLR\tbF\tbR\n')
% for i = 1:length(x)
%     fprintf('%.4f\t%.4f\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t%.4f\t%.4f\t%.4f \n',x(i),KxVec(i),tan(phiF(i)),tan(phiR(i)),sigmaCLF(i),sigmaCLR(i),bCLF(i),bCLR(i), bF(i), bR(i))
% end
%%
out.pitch = pitch;
out.pitch75 = pitch75;
out.xi = x;
out.cF = bF;
out.cR = bR;
out.betaF = betaF;
out.betaR = betaR;
out.betadegF = rad2deg(betaF);
out.betadegR = rad2deg(betaR);
out.beta75 = beta75;
out.betadeg75 = betadeg75;
out.beta75F = beta75F;
out.beta75R = beta75R;
out.sigmaF = sigmaF;
out.sigmaR = sigmaR;
out.sigmaCLF = sigmaCLF;
out.sigmaCLR = sigmaCLR;
out.CLF = CLF;
out.CLR = CLR;
out.cCLF = bCLF;
out.cCLR = bCLR;
% out.CtF = CtF;
% out.CtR = CtR;
out.Ct = Ct;
% out.CpF = CqF*omega;
% out.CpR = CqR*omega;
out.Cp = Cp;
% out.CqF = CqF;
% out.CqR = CqR;
out.Cq = Cq;
% out.TF = TF;
% out.TR = TR;
out.T = Ts;
% out.Tgil = Tgil;
% out.PF = PF;
% out.PR = PR;
out.P = P;
% out.QF = QF;
% out.QR = QR;
out.Q = Q;
out.J = J;
% out.etaF = etaF;
% out.etaR = etaR;
out.eta = eta_D;
% out.etagil = etagil;
out.TAFF = TAFF;
out.TAFR = TAFR;
out.TAF = TAF;
out.solidityF = solidityF;
out.solidityR = solidityR;
out.solidity = solidity;
