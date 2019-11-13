
% Wdot = 900*550;
% rho = 1.225;
% V = 625;
% n = 25;
% D = 7.5;
% B = 4;
% h = 28000*.305;
% 
% Wdot = 400e3;
% rho = 1.225;
% V = 153;
% n = 47;
% D = 1.53;
% B = 4;
% h = 0;

% Liebeck Single
clear;
Wdot = 70*550;
rho = 0.00237717;
V = 161.33;
h=0;
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
CL_design = 0.7;
xi = [0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 0.99];

design_ls = liebeckSingle(inputs,polares,propType,CL_design,xi)
% lscs = liebeckSingle_perf(inputs,design_cs)
lsls = liebeckSingle_perf(inputs,design_ls)


% des = davidsonDual(inputs,polares,3,0.52,xi)
% % an = davidsonDual_perf(inputs,des)
% 
% 
% figure(1)
% hold on
% plot(des.xi,des.cF,des.xi,des.cR)
% xlabel('x = r/R')
% ylabel('Chord [m]')
% %%
% figure(2)
% hold on
% plot(des.xi,des.betadeg)
% xlabel('x = r/R')
% ylabel('\beta [°]')
% 
% figure(3)
% hold on
% plot(an.J,an.Ct2)
% xlabel('J')
% ylabel('C_t')
% 
% figure(4)
% hold on
% plot(an.J,an.Cp2)
% xlabel('J')
% ylabel('C_p')
% 
% figure(5)
% hold on
% plot(an.J,an.eta)
% xlabel('J')
% ylabel('\eta')
% 
% figure(6)
% hold on
% plot(an.J,an.T)
% xlabel('J')
% ylabel('T')


