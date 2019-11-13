clear; clc; clear global
% AZX02
% global Wdot rho n D B R omega rpm
Wdot = 200e3;
h = 0;
V = 153;
n = 47;
D = 1.53;
B = 4;

crieglerPlots;

[T, a, patm, rho] = atmosisa(h);

rpm = 60*n;
J = V/(n*D);
nP = 'dualRotation'; % 'singleRotation' or 'dualRotation'
interpMode = 'spline'; % 'linear' or 'spline
R = D/2;
omega = n*2*pi;
F = pi*R^2; % Projected area of the helix
W = sqrt(V^2 + (omega*R)^2);

Cp = Wdot/(rho*n^3*D^5);
Cq = Wdot/(2*pi*rho*n^2*D^5);
% Pct = Wdot/(0.5*rho*V^3*pi*(D/2)^2);

load('propGeomCrigler.mat')

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
    load('polares360.mat');
    pol.CL_CD = pol.CL./pol.CD;
end

%%

V = 5:5:180;
% V = 130;
for jj = 1:length(V)
    J(jj) = V(jj)/(n*D);
    
    vr = sqrt(V(jj)^2 + (omega.*x.*R).^2);
    nu = 0.000014607;
    Re = vr.*bF/nu;
    Mach = vr/340;
    
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
    phiF = atan(J(jj)./(pi.*x).*(1+0.5*wbar*(1+0.5*k.*(tan(phi).^2))));
    phiR = atan(J(jj)./(pi.*x).*(1+0.5*wbar*(1-0.5*k.*(tan(phi).^2))));
    x1 = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
    [KxVec1] = getKx(nP, interpMode, B, Jw(jj)); % KxVec comes for x1 = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
    KxVec = interp1(x1,KxVec1,x,'linear', 'extrap');
    
    erro = inf;
    %%
    while erro > 1e-4
        
        alphaF = deg2rad(betaF') - phiF;
        alphaR = deg2rad(betaR') - phiR;
        
        for i = 1:length(x)
            CLFnew(i) = interp1(pol.alpha,pol.CL,rad2deg(alphaF(i)));
            CLRnew(i) = interp1(pol.alpha,pol.CL,rad2deg(alphaR(i)));
            CDF(i) = interp1(pol.alpha,pol.CD,rad2deg(alphaF(i)));
            CDR(i) = interp1(pol.alpha,pol.CD,rad2deg(alphaR(i)));
        end
        
        erro = abs(CLFnew(i) - CLF(i))
        
        CLF = CLFnew;
        CLR = CLRnew;
    end
    %%
    
    sigmaCLF = sigmaF.*CLF;
    sigmaCLR = sigmaR.*CLR;
    
    gammaF = atan(CDF./CLF);
    gammaR = atan(CDR./CLR);
    dCqdxF = pi/8*J(jj)^2*(1+1/4*k*wbar*sin(phi0).^2).^2.*sigmaCLF.*cos(phiF)./(sin(phi0).^2).*(tan(phiF) + tan(gammaF)).*x.^2;
    CqF = trapz(x,dCqdxF);
    dCqdxR = pi/8*J(jj)^2*(1+3/4*k*wbar*sin(phi0).^2).^2.*sigmaCLR.*cos(phiR)./(sin(phi0).^2).*(tan(phiR) + tan(gammaR)).*x.^2;
    CqR = trapz(x,dCqdxR);
    dCtdxF = 2./x.*(1 - tan(phiF).*tan(gammaF))./(tan(phiF) + tan(gammaF)).*dCqdxF;
    CtF = trapz(x,dCtdxF);
    dCtdxR = 2./x.*(1 - tan(phiR).*tan(gammaR))./(tan(phiR) + tan(gammaR)).*dCqdxR;
    CtR = trapz(x,dCtdxR);
    Cq = CqF + CqR;
    Cp(jj) = 2*pi*Cq;
    Ct(jj) = CtF + CtR;
    
    eta(jj) = J(jj)*Ct(jj)/Cp(jj)
    T(jj) = Ct(jj)*(rho*n^2*D^4)
    
    
    %%
end
%%

hold on
figure(1)
subplot(2,2,1)
hold on
plot(J,eta)
xlabel('J')
ylabel('\eta')
subplot(2,2,2)
hold on
plot(J,Cp)
xlabel('J')
ylabel('C_p')
subplot(2,2,3)
hold on
plot(J,Ct)
xlabel('J')
ylabel('C_t')
subplot(2,2,4)
hold on
plot(J,T)
xlabel('J')
ylabel('T [N]')




