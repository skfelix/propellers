function analysis = liebeckSingle_perf(inputs,design)

Wdot = inputs.Wdot;
V = inputs.V;
n = inputs.n;
D = inputs.D;
rho = inputs.rho;
B = inputs.B;
h = inputs.h;

xi = design.xi;
c  = design.c;
beta(1,:) = design.beta;
sigma = design.sigma;

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
% V=153;
% V = (45:5:185)/.305;
r = xi*R;

for i = 1:length(V)
    
    J(i) = V(i)/(n*D);
    omega = n*2*pi;
    lambda = V(i)/(omega*R);
    x = omega.*r/V(i);
    Pc1 = Wdot/(0.5*rho*V(i)^3*pi*R^2);
    Cp1 = Wdot/(rho*n^3*D^5);    
    
    phi = atan(V(i)./(omega*r));
    erro = inf;
    
    while erro > 1e-2
        
        phit = atan(tan(phi).*xi);
        f = 0.5*B*(1-xi)./sin(phit);
        F = (2/pi*acos(exp(-f)));
        Fh = 2/pi*acos(exp(-B/2*(x-0.2)./(0.2*sin(phi))));
        F = F.*Fh;
        G = F.*x.*cos(phi).*sin(phi);
      
        alpha = beta - phi;
        
        for j = 1:length(r)
            CL(j) = interp1(pol.alpha,pol.CL,rad2deg(alpha(j)));
            CD(j) = interp1(pol.alpha,pol.CD,rad2deg(alpha(j)));
        end
        
        e = CD./CL;
        Cy = CL.*(cos(phi) - e.*sin(phi));
        Cx = CL.*(sin(phi) + e.*cos(phi));
        
        Ka = Cy./(4*sin(phi).^2);
        Kb = Cx./(4*cos(phi).*sin(phi));        
        
        a = sigma.*Ka./(F - sigma.*Ka);
        b = sigma.*Kb./(F + sigma.*Kb);
        a(end) = 0;
        b(end) = 0;

        phinew = atan(V(i).*(1+a)./(omega.*r.*(1-b)));
       
        erro = norm(abs(phi-phinew))
        
        k = 0.9;
        phi = (1-k)*phi + k*phinew;        
    end   
   
    W = V(i).*(1+a)./sin(phi);
    
    dT = 0.5*rho.*W.^2*B.*c.*Cy;
    T(i) = trapz(r,dT)
    dQ = 0.5*rho.*W.^2*B.*c.*Cx.*r;
    Q = trapz(r,dQ);
    P = omega*Q;    
    Ct(i) = T(i)./(rho*n^2*D^4);
    Cp(i) = P./(rho*n^3*D^5);
    eta(i) = J(i)*Ct(i)/Cp(i); 
    
%     dCt = pi^3/4*sigma.*Cy.*xi.*F.^32./((F+sigma.*Kb).*cos(phi)).^2;
%     dCt(end) = 0;
    dCt = pi^3/4*sigma.*Cy.*xi.^3.*F.^2./((F+sigma.*Kb).*cos(phi)).^2;    
    Ct2(i) = real(trapz(xi,dCt));
    dCp = pi*dCt.*xi.*Cx./Cy;
    Cp2(i) = real(trapz(xi,dCp));
end

analysis.J = J;
analysis.Ct = Ct;
analysis.Cp = Cp;
analysis.Ct2 = Ct2;
analysis.Cp2 = Cp2;
analysis.eta = eta;
analysis.T = T;