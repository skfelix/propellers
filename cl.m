% clear all global; clc
% AZX02
% clear;
Wdot = 400e3;
rho = 1.225;
V = 153;
n = 47;
D = 1.53;
B = 4;
h = 0;

% Crigler
% clear;
% Wdot = 2000*550;
% rho = 0.001065;
% V = 623;
% n = 23;
% D = 12;
% B = 4;
% h = 0;

% Wdot = 2000*550;
% rho = 0.001065;
% V = 623;
% n = 23;
% D = 12;
% B = 4;
% h=0;


% Liebeck Single
% clear;
% Wdot = 70*550;
% rho = 0.00237717;
% V = 161.33;
% n = 40;
% D = 5.75;
% B = 2;

% Liebeck Single
% clear
% Wdot = 70*745;
% Wdot = 53e3;
% rho = 1.225;
% V = 49.17;
% n = 40;
% D = 1.753;
% B = 2;

% Liebeck Dual
% clear
% Wdot = 150*550;
% rho = 0.00237717;
% V = 225;
% n = 45;
% D = 5.75;
% B = 6;
% h=0;

inputs.Wdot = Wdot;
inputs.V = V;
inputs.rho = rho;
inputs.n = n;
inputs.D = D;
inputs.B = B;
inputs.h = h;

polares =0;
propType = 3;
xi = [0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 0.99];


des1 = liebeckSingle(inputs,polares,1,1,xi)
des2 = liebeckSingle(inputs,polares,2,1,xi)
des3 = liebeckSingle(inputs,polares,3,0.5,xi)
des4 = liebeckSingle(inputs,polares,3,1.2,xi)

an1 = liebeckSingle_perf(inputs,des1)
an2 = liebeckSingle_perf(inputs,des2)
an3 = liebeckSingle_perf(inputs,des3)
an4 = liebeckSingle_perf(inputs,des4)

figure(1)
subplot(1,2,1)
title('a)','FontWeight','bold')
hold on
plot(des2.xi,des2.c,'ok-','LineWidth',1.2)
plot(des3.xi,des3.c,'^k-','LineWidth',1.2)
plot(des1.xi,des1.c,'sk-','LineWidth',1.2)
plot(des4.xi,des4.c,'dk-','LineWidth',1.2)
xlabel('x = r/R','FontWeight','bold')
ylabel('c [m]','FontWeight','bold')
legend('C_L = 0,4116 (C_D)_{min}','CL = 0,5','C_L = 0,7953 (C_L/C_D)_{max}', 'C_L = 1,2','Location','best','FontWeight','bold')
grid on; grid minor;

subplot(1,2,2)
title('b)','FontWeight','bold')
hold on
plot(des1.xi,des1.betadeg,'sk-','LineWidth',1.2)
plot(des2.xi,des2.betadeg,'ok-','LineWidth',1.2)
plot(des3.xi,des3.betadeg,'^k-','LineWidth',1.2)
plot(des4.xi,des4.betadeg,'dk-','LineWidth',1.2)
xlabel('x = r/R','FontWeight','bold')
ylabel('\beta [°]','FontWeight','bold')
legend('C_L = 0,4116 (C_D)_{min}','CL = 0,5','C_L = 0,7953 (C_L/C_D)_{max}', 'C_L = 1,2','Location','best','FontWeight','bold')
grid on; grid minor;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)

subplot(2,2,1)
title('a)','FontWeight','bold')
hold on
plot(an2.J,an2.eta,'ok-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an3.J,an3.eta,'^k-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an1.J,an1.eta,'sk-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an4.J,an4.eta,'dk-','LineWidth',1.2,'MarkerIndices',1:3:27)
xlabel('J = V/nD','FontWeight','bold')
ylabel('\eta','FontWeight','bold')
legend('C_L = 0,4116 (C_D)_{min}','CL = 0,5','C_L = 0,7953 (C_L/C_D)_{max}', 'C_L = 1,2','Location','best','FontWeight','bold')
grid on; grid minor;

subplot(2,2,2)
title('b)','FontWeight','bold')
hold on
plot(an2.J,an2.Cp,'ok-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an3.J,an3.Cp,'^k-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an1.J,an1.Cp,'sk-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an4.J,an4.Cp,'dk-','LineWidth',1.2,'MarkerIndices',1:3:27)
xlabel('J = V/nD','FontWeight','bold')
ylabel('C_P','FontWeight','bold')
legend('C_L = 0,4116 (C_D)_{min}','CL = 0,5','C_L = 0,7953 (C_L/C_D)_{max}', 'C_L = 1,2','Location','best','FontWeight','bold')
grid on; grid minor;

subplot(2,2,3)
title('c)','FontWeight','bold')
hold on
plot(an2.J,an2.Ct,'ok-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an3.J,an3.Ct,'^k-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an1.J,an1.Ct,'sk-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an4.J,an4.Ct,'dk-','LineWidth',1.2,'MarkerIndices',1:3:27)
xlabel('J = V/nD','FontWeight','bold')
ylabel('C_T','FontWeight','bold')
legend('C_L = 0,4116 (C_D)_{min}','CL = 0,5','C_L = 0,7953 (C_L/C_D)_{max}', 'C_L = 1,2','Location','best','FontWeight','bold')
grid on; grid minor;

subplot(2,2,4)
title('d)','FontWeight','bold')
hold on
plot(an2.J,an2.T,'ok-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an3.J,an3.T,'^k-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an1.J,an1.T,'sk-','LineWidth',1.2,'MarkerIndices',1:3:27)
plot(an4.J,an4.T,'dk-','LineWidth',1.2,'MarkerIndices',1:3:27)
xlabel('J = V/nD','FontWeight','bold')
ylabel('T [N]','FontWeight','bold')
legend('C_L = 0,4116 (C_D)_{min}','CL = 0,5','C_L = 0,7953 (C_L/C_D)_{max}', 'C_L = 1,2','Location','best','FontWeight','bold')
grid on; grid minor;