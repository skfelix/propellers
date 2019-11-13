function analysis = crieglerSingle_perf(inputs,design)

Wdot = inputs.Wdot;
V = inputs.V;
n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;

x = design.xi;
b  = design.c;
beta(1,:) = design.beta;
sigma = design.sigma;
CL = design.CL;

crieglerPlots;

[T, asound, patm, rho] = atmosisa(h);
R = D/2;
nu = 0.000014607;
rpm = 60*n;
omega = n*2*pi;
J = V/(n*D);
lambda = V/(omega*R);
S = pi*R^2; % Projected area of the helix
nP = 'singleRotation'; % 'singleRotation' or 'dualRotation'
interpMode = 'spline'; % 'linear' or 'spline

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

%%

V = 50:5:180;
% V = 153;
for jj = 1:length(V)
    J(jj) = V(jj)/(n*D);
    
    W = sqrt(V(jj)^2 + (omega.*x.*R).^2);
    nu = 0.000014607;
    Re = W.*b/nu;
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
    
    erro = inf;
    %%
    while erro > 1e-4
        
        alpha = beta - phi;
        
        for i = 1:length(x)
            CLnew(i) = interp1(pol.alpha,pol.CL,rad2deg(alpha(i)));
            CD(i) = interp1(pol.alpha,pol.CD,rad2deg(alpha(i)));
        end
        
        erro = abs(CLnew - CL)        
        CL = CLnew;
    end

%     p = 0.5.*sin(phi0).*cos(phi0);
%     q = 0.5./sin(phi0);
%     r = cos(phi0).^2 - sin(phi0).^2;
%     s = Jw(jj)./(pi*x).*(KxVec.*sin(phi0));
%     
%     
%     X0 = q.*s./(2.*p - q.*r.*s);
%     b_line = 1./(4*X0.*sin(phi0));
%     
%     while erro > 1e-2
%         
%         alphai = b_line.*sigma.*CL.*(1+X0.*cos(phi0).^2);
%         alpha = beta - phi0 - alphai;
%         
%         for i = 1:length(x)
%             CLnew(i) = interp1(pol.alpha,pol.CL,rad2deg(alpha(i)));
%             CD(i) = interp1(pol.alpha,pol.CD,rad2deg(alpha(i)));
%         end
%         
%         erro = abs(CLnew(i) - CL(i))
%         
%         CL = CLnew;
%     end
    %%
    
    sigmaCL = sigma.*CL;

 
    e = atan(CD./CL);
    dCqdx = pi/8*J(jj)^2*((1+0.5*wbar*cos(phi).^2)./sin(phi)).^2.*sigmaCL.*(sin(phi) + e.*cos(phi)).*x.^2;
    
    Cq = trapz(x,dCqdx);
    % dCtdx = 2./x.*(1 - tan(phi).*tan(e))./(tan(phi) + tan(e)).*dCqdx;
    dCtdx = pi/4*J(jj)^2*((1+0.5*wbar*cos(phi).^2)./sin(phi)).^2.*sigmaCL.*(cos(phi) - e.*sin(phi)).*x;
    Ct(jj) = trapz(x,dCtdx);
    Cp(jj) = 2*pi*Cq
    
    eta(jj) = J(jj)*Ct(jj)/Cp(jj);
    T(jj) = Ct(jj)*(rho*n^2*D^4);
    
end


%%
analysis.J = J;
analysis.Ct = Ct;
analysis.Cp = Cp;
analysis.eta = eta;
analysis.T = T;



