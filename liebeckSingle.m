% function out = liebeckSingle(inputs,polares,propType,CL_design,xi)

h=0;
Wdot = 70*550;
rho = 0.00237717;
V = 161.33;
n = 40;
D = 5.75;
B = 2;

inputs.Wdot = Wdot;
inputs.V = V;
inputs.rho = rho;
inputs.n = n;
inputs.D = D;
inputs.B = B;
inputs.h = h;

polares =0;
propType = 3;
% CL_design = 0.5;
CL_design = 0.7;
xi = [0.173913043	0.311582609	0.449286957	0.586956522	0.724626087	0.862330435	1];

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
Pc = Wdot/(0.5*rho*V^3*pi*R^2);

%% Polares - Executar somente se mudar as condições de contorno
% polares = 0; % 0 para não gerar polar, 1 para gerar
if (polares == 1)
    [p f] = xfoil('clark-y.dat',0:0.2:8,0.5e6,0,'oper iter 100');
    p.CL_CD = p.CL./p.CD;
%     p.alpha = unique(p.alpha,'stable'
    save('polar_clarky.mat','p');
else
    load('polar_clarky.mat');
end


%% Prandtl Factor

% xi = [0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 0.99];
% xi = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99];
% % xi = [0.5 0.8958 1.2917 1.6875 2.0833 2.4792 2.875]*.3048/R;
% xi = linspace(0.5,2.875,21)/R;
% xi = 0.1:0.02:0.99;
r = xi*R;
W = sqrt(V^2+(omega*r).^2);
% Mach = W/asound;


zeta = 0.2;
erro = inf;

while erro > 1e-4 
    
    phit = atan(lambda*(1+0.5*zeta));
    
    x = omega*r/V;
    f = 0.5*B*(1-xi)/sin(phit);
    F = 2/pi*acos(exp(-f));
    phi = atan(tan(phit)./xi);
    Fh = 2/pi*acos(exp(-B/2*(x-0.2)./(0.2*sin(phi))));
    F = F.*Fh;
    G = F.*x.*cos(phi).*sin(phi);

    
    %% Prop Type
    
    if propType == 1 % max CL/CD
        for i = 1:length(r)
            idx = find(p.CL_CD==max(p.CL_CD));
            alpha(i) = p.alpha(idx);
            CL(i) = p.CL(idx);
            CD(i) = p.CD(idx);
        end
        e = 1/p.CL_CD(idx);
    elseif propType == 2 % min CD
        for i = 1:length(r)
            idx = find(p.CD==min(p.CD));
            alpha(i) = p.alpha(idx(1));
            CL(i) = p.CL(1);
            CD(i) = p.CD(1);
        end
        e = CD./CL;
    elseif propType == 3 % fixed CL | CL design
        for i = 1:length(r)
            alpha(i) = interp1(p.CL,p.alpha,CL_design,'linear','extrap');
            CL(i) = CL_design;
            CD(i) = interp1(p.alpha,p.CD,alpha(i),'linear','extrap');
        end
        e = CD./CL;
    elseif propType == 4 % Takeoff
        for i = 1:length(r)
            alpha(i) = interp1(p.CL,p.alpha,CL_design,'linear','extrap');
            CL(i) = CL_design;
            CD(i) = interp1(p.alpha,p.CD,alpha(i),'linear','extrap');
        end
        e = CD./CL;
    end
        
    Cy = CL.*(cos(phi) - e.*sin(phi));
    Cx = CL.*(sin(phi) + e.*cos(phi));
    
    %%
    
    Wc = 4*pi*lambda*G*V*R*zeta./(CL*B);
    Re = Wc/nu;
    
    a = 0.5*zeta*cos(phi).^2.*(1-e.*tan(phi));
    b = 0.5*zeta.*x.*cos(phi).*sin(phi).*(1+e./tan(phi)); 
%     a(end) = 0;
%     b(end) = 0;
    W = V*(1+a)./sin(phi);
    c = Wc./W;
    beta = deg2rad(alpha)+phi;
    
    dI1 = 4*xi.*G.*(1-e.*tan(phi));
    dI2 = lambda.*(0.5*dI1./xi).*(1+e./tan(phi)).*sin(phi).*cos(phi);
    dJ1 = 4*xi.*G.*(1+e./tan(phi));
    dJ2 = 0.5*dJ1.*(1-e.*tan(phi)).*cos(phi).^2;
    
    I1 = trapz(xi,dI1);
    I2 = trapz(xi,dI2);
    J1 = trapz(xi,dJ1);
    J2 = trapz(xi,dJ2);
    zetanew = -(0.5*J1./J2) + sqrt((0.5*J1./J2).^2 + Pc./J2);
    
    Tc = I1.*zeta - I2.*zeta.^2;
   
    erro = abs(zeta - zetanew);
    
    zeta = zetanew;
    
end

Jw = mean(V*(1+a)/n/D);
sigma = B*c./(2*pi.*r);
sigmaCL = sigma.*CL;
cCL = c.*CL;

% 
eta = Tc/Pc;
dT = 0.5*rho.*W.^2*B.*c.*Cy;
T = trapz(r,dT);
dQ = 0.5*rho.*W.^2*B.*c.*Cx.*r;
Q = trapz(r,dQ);
P = omega*Q;
Ct = T/(rho*n^2*D^4);
Cp = P/(rho*n^3*D^5);
betadeg = rad2deg(beta);
phideg = rad2deg(phi);
beta75 = interp1(xi,beta,0.75);
betadeg75 = rad2deg(beta75);
sigma75 = B*interp1(xi,c,0.75)/(2*pi*0.75*R);
AF = 1e5/(32*R^5)*trapz(xi*R, c.*(xi*R).^3);
TAF = B*AF;
pitch = 2*pi*(xi*R).*tan(beta);
pitch75 = 2*pi*(0.75*R)*tan(beta75);

%%

deltapa = trapz(xi*D,0.5*rho*((V*(1+a)).^2-V^2))
deltapr = trapz(r,rho*(omega-0.5*b*omega).*b.*omega.*r.^2)% + ...
%           trapz(r,0.5*rho*(omega^2*r.^2))
deltap = deltapa - deltapr
Ta = trapz(r,4*rho*pi*r.*V^2.*a.*(1+a).*F)
trapz(r,rho*pi*r.*V^2.*a.*(2+3*a).*F) + deltap*pi*R^2

%%

out.pitch = pitch;
out.pitch75 = pitch75;
out.xi = xi;
out.c = c;
out.cCL = cCL;
out.beta = beta;
out.betadeg = betadeg;
out.beta75 = beta75;
out.betadeg75 = betadeg75;
out.sigma = sigma;
out.sigmaCL = sigmaCL;
out.cCL = cCL;
out.Ct = Ct;
out.Cp = Cp;
out.CL = CL;
out.Cq = Cp/omega;
out.T = T;
out.P = P;
out.Q = Q;
out.J = J;
out.eta = eta;
out.sigma75 = sigma75;
out.TAF = TAF;

