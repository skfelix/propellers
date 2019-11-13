function out = crieglerSingle(inputs,polares,propType,CL_design,xi)

Wdot = inputs.Wdot;
V = inputs.V;
n = inputs.n;
rho = inputs.rho;
D = inputs.D;
B = inputs.B;
h = inputs.h;

crieglerPlots;

% [T, asound, patm, rho] = atmosisa(h);
nu = 0.000014607;
rpm = 60*n;
J = V/(n*D);
nP = 'singleRotation'; % 'singleRotation' or 'dualRotation'
interpMode = 'spline'; % 'linear' or 'spline
R = D/2;
J = V/(n*D);
omega = n*2*pi;
F = pi*R^2; % Projected area of the helix
W = sqrt(V^2 + (omega*R)^2);

Cpt = Wdot/(rho*n^3*D^5);
Pct = Wdot/(0.5*rho*V^3*pi*R^2);
%%

for wbar = 0:0.00001:1
    Jw = J*(1+wbar);
    [k,ek] = getKappa(nP, interpMode, B, Jw);
    Pc = 2*k*wbar*(1+wbar)*(1+ek*wbar);
    if abs(0.988*Pct - Pc) < 1e-4
        break;
    end
end

Pc_k = Pc/k;
eta = getEff(Pc_k, ek, interpMode);

%% Tables
x = xi;
w = wbar*V;
x_Kx = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
r = x*R;
[KxVec1] = getKx(nP, interpMode, B, Jw); % KxVec comes for x1 = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
% if x ~ the above, it is necessary to interpolate
KxVec = interp1(x_Kx,KxVec1,x,'linear', 'extrap');
phi = atan(1/pi*J*(1+0.5*wbar)./x);
phi0 = atan(J./(pi.*x));
sigmaCL = (1+wbar)./((1+0.5*wbar).*(1+0.5.*wbar.*cos(phi).^2)).*2*wbar.*KxVec.*sin(phi).^2./cos(phi);
bCL = sigmaCL.*2.*pi.*(x.*R)./B;

CL_design = 0.5;
b_guess = bCL/CL_design;
nu = 0.000014607;
vr = sqrt(V^2 + (omega.*x.*R).^2);
Re = vr.*b_guess/nu;
Mach = vr/340;

%% Polares - Executar somente se mudar as condições de contorno
% polares =0 ; % 0 para não gerar polar, 1 para gerar
if (polares == 1)
    %     [pol foil] = xfoil('clark-y.dat',-5:0.1:15,mean(Re),0.3,'oper iter 300');
    for i = 1:length(x)
        [p f] = xfoil('clark-y.dat',0:0.2:8,5e5,0.3,'oper iter 100');
        p.CL_CD = p.CL./p.CD;
        polRe{i} = p;
        foilRe{i} = f;
    end
    %     save('polares.mat','polRe','foilRe','pol','foil');
    save('polar_clarky.mat','polRe','foilRe');
else
    load('polar_clarky.mat');
end

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
    for i = 1:length(r)
        alpha(i) = interp1(p.CL,p.alpha,CL_design);
        CL(i) = CL_design;
        CD(i) = interp1(p.alpha,p.CD,alpha(i));
    end
    e = CD./CL;
elseif propType == 4 % Takeoff
    for i = 1:length(r)
        alpha(i) = interp1(p.CL,p.alpha,CL_design);
        CL(i) = CL_design;
        CD(i) = interp1(p.alpha,p.CD,alpha(i));
    end
    e = CD./CL;
end

%%
betadeg = rad2deg(phi)' + alpha';
beta = deg2rad(betadeg');
beta75 = interp1(x,beta,0.75);
betadeg75 = rad2deg(beta75);
sigma = sigmaCL./CL;
b = bCL./CL;


pitch = 2*pi*(x*R).*tan(beta);
pitch75 = 2*pi*(0.75*R)*tan(beta75);


lambda_g = J/pi;
sigmaCD = sigma.*CD;
tr = 2/lambda_g^2*trapz(x,sigmaCD.*x.^3./sin(phi)); % rotational-drag-loss coeff (torque loss)
ta = 2*trapz(x,sigmaCD.*x./sin(phi)); % axial-drag-loss coeff (thrust loss)
cs = 2*k*wbar*(1 + wbar*(0.5 + ek)); % total thrust coeff
eta_i = cs/Pct; % no drag efficiency
csT = cs - ta; % net thrust coeff
Pci = Pc + tr;
eta_D = csT/Pci;
T = rho*F*k*w*(V+w*(0.5+ek)); % total thrust
Ts = 1/2*rho*V^2*F*csT; % net thrust
sigma75 = B*interp1(x,b,0.75)/(2*pi*0.75*R);
AF = 1e5/(32*R^5)*trapz(x*R, b.*(x*R).^3);
TAF = B*AF;
% Ct = T/(rho*n^2*D^4);
% Cq = Cpt/(2*pi);
% Q = Cq*rho*n^2*D^5;

e = atan(CD./CL);
dCqdx = pi/8*J^2*((1+0.5*wbar*cos(phi).^2)./sin(phi)).^2.*sigmaCL.*(sin(phi) + e.*cos(phi)).*x.^2;

Cq = trapz(x,dCqdx);
% dCtdx = 2./x.*(1 - tan(phi).*tan(e))./(tan(phi) + tan(e)).*dCqdx;
dCtdx = pi/4*J^2*((1+0.5*wbar*cos(phi).^2)./sin(phi)).^2.*sigmaCL.*(cos(phi) - e.*sin(phi)).*x;
Ct = trapz(x,dCtdx);
Cp = 2*pi*Cq;

eta = J*Ct/Cp
T = Ct*(rho*n^2*D^4)
%%
% fprintf('=================== Propeller Characteristics ===================\n')
%     fprintf('x\t\ttan(phi)\tK(x)\tsigmaCL\tbCL\t\tb\t\tb_falk\n')
%     for i = 1:length(x)
%         fprintf('%.4f\t%.4f\t\t%.4f\t%.4f\t%.4f\t%.4f \n',x(i),tan(phi(i)),KxVec(i), sigmaCL(i), bCL(i), b(i))
%     end

%%

out.pitch = pitch;
out.pitch75 = pitch75;
out.xi = x;
out.c = b;
out.beta = beta;
out.betadeg = betadeg;
out.beta75 = beta75;
out.betadeg75 = betadeg75;
out.sigma = sigma;
out.sigmaCL = sigmaCL;
out.CL = CL;
out.cCL = bCL;
out.Ct = Ct;
out.Cp = Cpt;
out.Cq = Cq;
out.T = T;
% out.P = P;
% out.Q = Q;
out.J = J;
out.eta = eta;
out.sigma75 = sigma75;
out.TAF = TAF;

