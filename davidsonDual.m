function out = davidsonDual(inputs,polares,propType,CL_design,xi)
% Wdot = 900*550;
% rho = 1.225;
% V = 625;
% n = 25;
% D = 7.5;
% B = 4;
% h = 28000*.305;
% 
% 
% inputs.Wdot = Wdot;
% inputs.V = V;
% inputs.rho = rho;
% inputs.n = n;
% inputs.D = D;
% inputs.B = B;
% inputs.h = h;
% 
% polares =0;
% propType = 3;
% CL_design = 0.52;
% xi = [0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 0.99];



Wdot = inputs.Wdot;
V = inputs.V;
n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;

crieglerPlots;

[T, asound, patm, rho] = atmosisa(h);
% rho = rho*0.00194032;
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

Cp = Wdot/(rho*n^3*D^5);
Pct = Wdot/(0.5*rho*V^3*pi*(D/2)^2);

%% Polares
% polares = 0; % 0 para não gerar polar, 1 para gerar
if (polares == 1)
    [pol f] = xfoil('clark-y.dat',0:0.05:8,0.5e6,0,'oper iter 100');
    pol.CL_CD = pol.CL./pol.CD;
    save('polar_clarky.mat','p');
else
    load('polar_clarky.mat');
end
pol = p;
%% Displacement Velocity w

% Crigler

for wbar = 0:0.0001:1
    Jw = V/(n*D)*(1+wbar);
    [masscoefk,ek] = getKappa(nP, interpMode, B, Jw);
    Pc = 2*masscoefk*wbar*(1+wbar)*(1+ek*wbar);
    if abs(Pct - Pc) < 1e-4
        break;
    end
end

% Davidson

% w = sqrt(2*patm/rho*((1+0.2*Minf^2)^3.5 -1));
% wbar = w/V;

w = wbar*V;
Jw = (V+w)/(n*D);
%% Displacement Velocity w

x1 = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
% x = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];
x = xi;
[KxVec1] = getKx(nP, interpMode, B, Jw); % KxVec comes for x1 = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
KxVec = interp1(x1,KxVec1,x,'linear', 'extrap');
phi = atan(1/pi*J*(1+0.5*wbar)./x);
phi0 = atan(J./(pi.*x));

Minf = V/asound;
vr = sqrt(V^2 + (omega.*x.*R).^2);
Mach = vr/340;

% Tip Loss Correction

p = 0.5.*sin(phi0).*cos(phi0);
q = 0.5./sin(phi0);
r = cos(phi0).^2 - sin(phi0).^2;
s = Jw./(pi*x).*(KxVec.*sin(phi0));


X0 = q.*s./(2.*p - q.*r.*s);

%% Prop Type
if propType == 1 % max CL/CD
    for i = 1:length(r)
        idx = find(pol.CL_CD==max(pol.CL_CD));
        alpha(i) = pol.alpha(idx);
        CL(i) = pol.CL(idx);
        CD(i) = pol.CD(idx);
    end
    e = 1/pol.CL_CD(idx);
elseif propType == 2 % min CD
    for i = 1:length(r)
        idx = find(pol.CD==min(pol.CD));
        alpha(i) = pol.alpha(idx(1));
        CL(i) = pol.CL(1);
        CD(i) = pol.CD(1);
    end
    e = CD./CL;
elseif propType == 3 % fixed CL | CL design
    for i = 1:length(r)
        alpha(i) = interp1(pol.CL,pol.alpha,CL_design,'linear','extrap');
        CL(i) = CL_design;
        CD(i) = interp1(pol.alpha,pol.CD,alpha(i),'linear','extrap');
    end
    e = CD./CL;
elseif propType == 4 % Takeoff
    for i = 1:length(r)
        alpha(i) = interp1(pol.CL,pol.alpha,CL_design,'linear','extrap');
        CL(i) = CL_design;
        CD(i) = interp1(pol.alpha,pol.CD,alpha(i),'linear','extrap');
    end
    e = CD./CL;
end

Cy = CL.*(cos(phi) - e.*sin(phi));
Cx = CL.*(sin(phi) + e.*cos(phi));

phiq0 = sin(phi0) + e.*cos(phi0);
%% Power coefficient == Propeller power absorted

for k = 0.1:0.00001:0.6
    sigmaCL = (0.5*k*phiq0.*cos(phi0) - e)./ ...
        (0.5./sin(phi0).*(1./X0 + cos(phi0).^2 - sin(phi0).^2));
    dCpdx = pi^4/4*sec(phi0).^2.*x.^4.*sigmaCL.*phiq0;
    Cp1 = trapz(x,dCpdx);
    if abs(Cp1 - Cp) < 1e-4
        break;
    end
end

b_line = 1./(4*X0.*sin(phi0));
A = pi*rho*(x*R).^4*omega^3.*sec(phi0).^3;
t1 = b_line.*(1 + X0.*cos(2*phi0));

etar = 1 - (t1.*(sigmaCL) + e).*sec(phi0)./(phiq0);
eta = mean(etar);


%% Geometry
sigmaCLF = sigmaCL/2;
sigmaCLR = sigmaCL/2;

bCLF = sigmaCLF.*2.*pi.*(x.*R)./(B/2);
bCLR = sigmaCLF.*2.*pi.*(x.*R)./(B/2);
sigmaF = sigmaCLF./CL;
sigmaR = sigmaCLR./CL;
bF = bCLF./CL;
bR = bCLR./CL;

alphaiF = b_line.*sigmaCLF.*(1+X0.*cos(phi0).^2);
alphaiR = b_line.*sigmaCLR.*(1+X0.*(cos(phi0).^2 - 2*sin(phi0).^2));
betaF = phi0 + deg2rad(alpha) + alphaiF;
betaR = phi0 + deg2rad(alpha) + alphaiR;
betadegF = rad2deg(betaF);
betadegR = rad2deg(betaR);
solidityF = B/2*interp1(x,bF,0.75)/(2*pi*0.75*R);
solidityR = B/2*interp1(x,bR,0.75)/(2*pi*0.75*R);
solidity = solidityF+solidityR;
TAFF = B/2*1e5/(32*R^5)*trapz(x*R, bF.*(x*R).^3);
TAFR = B/2*1e5/(32*R^5)*trapz(x*R, bR.*(x*R).^3);
TAF = TAFF + TAFR;
beta75F = interp1(x*R,betadegF,0.75);
beta75R = interp1(x*R,betadegR,0.75);



% eta = 1 - k/4;
Ct = eta*Cp/J;
T = Ct*(rho*n^2*D^4);


out.xi = x;
out.cF = bF;
out.cR = bR;
out.betaF = betaF;
out.betaR = betaR;
out.betadegF = rad2deg(betaF);
out.betadegR = rad2deg(betaR);
out.beta75F = beta75F;
out.beta75R = beta75R;
out.sigmaF = sigmaF;
out.sigmaR = sigmaR;
out.sigmaCLF = sigmaCLF;
out.sigmaCLR = sigmaCLR;
out.cCLF = bCLF;
out.cCLR = bCLR;
out.CLF = CL;
out.CLR = CL;
% out.CtF = CtF;
% out.CtR = CtR;
out.Ct = Ct;
% out.CpF = CpF;
% out.CpR = CpR;
out.Cp = Cp;
% out.CqF = CpF/omega;
% out.CqR = CpR/omega;
% out.Cq = Cp/omega;
% out.TF = TF;
% out.TR = TR;
out.T = T;
% out.PF = PF;
% out.PR = PR;
% out.P = P;
% out.QF = QF;
% out.QR = QR;
% out.Q = Q;
out.J = J;
% out.etaF = etaF;
% out.etaR = etaR;
out.eta = eta;
out.TAFF = TAFF;
out.TAFR = TAFR;
out.TAF = TAF;
out.solidityF = solidityF;
out.solidityR = solidityR;
out.solidity = solidity;


