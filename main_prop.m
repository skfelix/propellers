% clear all global; clc
% AZX02
% clear;
% Wdot = 400e3;
% rho = 1.225;
% V = 153;
% n = 47;
% D = 1.53;
% B = 4;
% h = 0;

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
% h=0;
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
clear
Wdot = 150*550;
rho = 0.00237717;
V = 225;
n = 45;
D = 5.75;
B = 6;
h=0;

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
% xi = [0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 0.99];
% xi = [0.173913043	0.311582609	0.449286957	0.586956522	0.724626087	0.862330435	1];
xi = [0.173913043	0.2152	0.256521739	0.297808696	0.339130435	0.380417391	0.42173913	0.463026087	0.504347826	0.545634783	0.586956522	0.628243478	0.669565217	0.710852174	0.752173913	0.79346087	0.834782609	0.876069565	0.917391304	0.958678261	1];

%% ############## ADICIONAR CL E XI ######################
roda =0;
if roda ==1
    design_ls = liebeckSingle(inputs,polares,propType,CL_design,xi)
    analysis_ls = liebeckSingle_perf(inputs,design_ls)
    
    design_ld = liebeckDual(inputs,polares,propType,CL_design,xi)
    analysis_ld = liebeckDual_perf(inputs,design_ld)
    
    design_cd = crieglerDual(inputs,polares,propType,CL_design,xi)
%     analysis_cd = crieglerDual_perf_Davidson(inputs,design_cd)
    analysis_cd = davidsonDual_perf(inputs,design_cd)
    
%     design_dd = davidsonDual(inputs,polares,propType,CL_design,xi)
    design_dd = playleDual(inputs,polares,propType,CL_design,xi)
    analysis_dd = davidsonDual_perf(inputs,design_dd)
    
    design_cs = crieglerSingle(inputs,polares,propType,CL_design,xi)
    analysis_cs = crieglerSingle_perf(inputs,design_cs)
    save('propdata','design_ls','analysis_ls','design_ld','analysis_ld','design_cs','analysis_cs','design_cd','analysis_cd','design_dd','analysis_dd')
else
    load('propdata')
end


%% eta, Ct, Cp, T
perf = 1;
if perf 
figure(1)
subplot(2,2,1)
title('a)','FontWeight','bold')
hold on
plot(analysis_cs.J,analysis_cs.eta,'k-o','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_ls.J,analysis_ls.eta,'r-s','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_cd.J,analysis_cd.eta,'b-d','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_dd.J,analysis_dd.eta,'c-^','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_ld.J,analysis_ld.eta,'m-v','LineWidth',1,'MarkerIndices',1:3:27)
xlabel('J = V/nD','FontWeight','bold')
ylabel('\eta','FontWeight','bold')
axis([0.5 2.6 0.4 1])
legend('Crigler Single','Liebeck Single','Crigler Dual','Davidson Dual', ...
    'Liebeck Dual','Location','best','FontWeight','bold')
grid minor; grid on
subplot(2,2,2)
title('b)','FontWeight','bold')
hold on
plot(analysis_cs.J,analysis_cs.Cp,'k-o','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_ls.J,analysis_ls.Cp,'r-s','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_cd.J,analysis_cd.Cp,'b-d','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_dd.J,analysis_dd.Cp,'c-^','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_ld.J,analysis_ld.Cp,'m-v','LineWidth',1,'MarkerIndices',1:3:27)
xlabel('J = V/nD','FontWeight','bold')
ylabel('C_P','FontWeight','bold')
axis([0.5 2.6 0 0.65])
legend('Crigler Single','Liebeck Single','Crigler Dual','Davidson Dual', ...
    'Liebeck Dual','Location','best','FontWeight','bold')
grid minor; grid on
subplot(2,2,3)
title('c)','FontWeight','bold')
hold on
plot(analysis_cs.J,analysis_cs.Ct,'k-o','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_ls.J,analysis_ls.Ct,'r-s','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_cd.J,analysis_cd.Ct,'b-d','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_dd.J,analysis_dd.Ct,'c-^','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_ld.J,analysis_ld.Ct,'m-v','LineWidth',1,'MarkerIndices',1:3:27)
xlabel('J = V/nD','FontWeight','bold')
ylabel('C_T','FontWeight','bold')
axis([0.5 2.6 0 0.4])
legend('Crigler Single','Liebeck Single','Crigler Dual','Davidson Dual', ...
    'Liebeck Dual','Location','best','FontWeight','bold')
grid minor; grid on
subplot(2,2,4)
title('d)','FontWeight','bold')
hold on
plot(analysis_cs.J,analysis_cs.T,'k-o','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_ls.J,analysis_ls.T,'r-s','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_cd.J,analysis_cd.T,'b-d','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_dd.J,analysis_dd.T,'c-^','LineWidth',1,'MarkerIndices',1:3:27)
plot(analysis_ld.J,analysis_ld.T,'m-v','LineWidth',1,'MarkerIndices',1:3:27)
xlabel('J = V/nD','FontWeight','bold')
ylabel('T [N]','FontWeight','bold')
axis([0.5 2.6 0 6e3])
legend('Crigler Single','Liebeck Single','Crigler Dual','Davidson Dual', ...
    'Liebeck Dual','Location','best','FontWeight','bold')
grid minor; grid on
end
%%
on = 1;
if on
figure(2)
subplot(3,2,1)
title('a)','FontWeight','bold')
hold on
plot(design_ls.xi,design_ls.c,'ro-','LineWidth',1)
plot(design_cs.xi,design_cs.c,'ks-','LineWidth',1)
xlabel('x = r/R','FontWeight','bold')
ylabel('Chord [m]','FontWeight','bold')
legend('Liebeck Single','Crigler Single','FontWeight','bold')
axis([design_ls.xi(1) 1 0 0.25])
grid minor; grid on

subplot(3,2,2)
title('b)','FontWeight','bold')
hold on
plot(design_cd.xi,design_cd.cF,'b','LineWidth',1)
plot(design_cd.xi,design_cd.cR,'-bd','LineWidth',1)
plot(design_dd.xi,design_dd.cF,'c-','LineWidth',1)
plot(design_dd.xi,design_dd.cR,'c-^','LineWidth',1)
plot(design_ld.xi,design_ld.cF,'m-','LineWidth',1)
plot(design_ld.xi,design_ld.cR,'m-v','LineWidth',1)
xlabel('x = r/R','FontWeight','bold')
ylabel('Chord [m]','FontWeight','bold')
lgd = legend('Crigler Dual F','Crigler Dual R','Davidson Dual F','Davidson Dual R','Liebeck Dual F','Liebeck Dual R','FontWeight','bold')
axis([design_ls.xi(1) 1 0 0.25])
grid minor; grid on

subplot(3,2,3)
title('c)','FontWeight','bold')
hold on
plot(design_ls.xi,design_ls.betadeg,'ro-','LineWidth',1)
plot(design_cs.xi,design_cs.betadeg,'ks-','LineWidth',1)
xlabel('x = r/R','FontWeight','bold')
ylabel('\beta [°]','FontWeight','bold')
legend('Liebeck Single','Crigler Single','FontWeight','bold')
grid minor; grid on

subplot(3,2,4)
title('d)','FontWeight','bold')
hold on
plot(design_cd.xi,design_cd.betadegF,'b','LineWidth',1)
plot(design_cd.xi,design_cd.betadegR,'-bd','LineWidth',1)
plot(design_dd.xi,design_dd.betadegF,'c-','LineWidth',1)
plot(design_dd.xi,design_dd.betadegR,'c-^','LineWidth',1)
plot(design_ld.xi,design_ld.betadegF,'m-','LineWidth',1)
plot(design_ld.xi,design_ld.betadegR,'m-v','LineWidth',1)
xlabel('x = r/R','FontWeight','bold')
ylabel('\beta [°]','FontWeight','bold')
lgd = legend('Crigler Dual F','Crigler Dual R','Davidson Dual F','Davidson Dual R','Liebeck Dual F','Liebeck Dual R','FontWeight','bold')
grid minor; grid on

subplot(3,2,5)
title('e)','FontWeight','bold')
hold on
plot(design_ls.xi,design_ls.sigmaCL,'ro-','LineWidth',1)
plot(design_cs.xi,design_cs.sigmaCL,'ks-','LineWidth',1)
xlabel('x = r/R','FontWeight','bold')
ylabel('\sigmaC_L [-]','FontWeight','bold')
legend('Liebeck Single','Crigler Single','FontWeight','bold')
grid minor; grid on
axis([design_ls.xi(1) 1 0 0.25])

subplot(3,2,6)
title('f)','FontWeight','bold')
hold on
plot(design_cd.xi,design_cd.sigmaCLF,'b','LineWidth',1)
plot(design_cd.xi,design_cd.sigmaCLR,'-bd','LineWidth',1)
plot(design_dd.xi,design_dd.sigmaCLF,'c-','LineWidth',1)
plot(design_dd.xi,design_dd.sigmaCLR,'c-^','LineWidth',1)
plot(design_ld.xi,design_ld.sigmaCLF,'m-','LineWidth',1)
plot(design_ld.xi,design_ld.sigmaCLR,'m-v','LineWidth',1,'MarkerIndices',1:3:27)
xlabel('x = r/R','FontWeight','bold')
ylabel('\sigmaC_L [-]','FontWeight','bold')
lgd = legend('Crigler Dual F','Crigler Dual R','Davidson Dual F','Davidson Dual R','Liebeck Dual F','Liebeck Dual R','FontWeight','bold')
axis([design_ls.xi(1) 1 0 0.25])
grid minor; grid on
end

%% Comparação entre métodos de análise

roda1 = 1;
if roda1
ddcd = davidsonDual_perf(inputs,design_cd)
ddld = davidsonDual_perf(inputs,design_ld)
dddd = davidsonDual_perf(inputs,design_dd)

ldcd = liebeckDual_perf(inputs,design_cd)
ldld = liebeckDual_perf(inputs,design_ld)
lddd = liebeckDual_perf(inputs,design_dd)

lscs = liebeckSingle_perf(inputs,design_cs)
lsls = liebeckSingle_perf(inputs,design_ls)

cscs = crieglerSingle_perf(inputs,design_cs)
csls = crieglerSingle_perf(inputs,design_ls)
save('propdata2','ddcd','ddld','dddd','ldcd','ldld','lddd','lscs','cscs','csls')
else
    load('propdata2')
end

%%
 figure
subplot(3,3,1)
title('a)','FontWeight','bold')
hold on
plot(lscs.J,lscs.eta,':o','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(lsls.J,lsls.eta,':s','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(cscs.J,cscs.eta,'-.o','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
plot(csls.J,csls.eta,'-.s','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
xlabel('J = V/nD','FontWeight','bold')
ylabel('\eta','FontWeight','bold')
lg = plot(NaN,NaN,'k-o',NaN,NaN,'k-s',NaN,NaN,':',NaN,NaN,'-.');
for i = 1:length(lg)
    lg(i).LineWidth = 1.2;
end
lg(3).Color = [0.9290 0.6940 0.1250];
lg(4).Color =[0.4940 0.1840 0.5560];
legend(lg,'Crigler Single','Liebeck Single','Liebeck Analysis','Davidson Mod. Analysis','Location','best','FontWeight','bold','Color','none')
axis([0.5 2.6 0.4 1])
grid minor; grid on

subplot(3,3,4)
title('d)','FontWeight','bold')
hold on
plot(ddcd.J,ddcd.eta,'-d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(dddd.J,dddd.eta,'-^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ddld.J,ddld.eta,'-v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ldcd.J,ldcd.eta,'--d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(lddd.J,lddd.eta,'--^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(ldld.J,ldld.eta,'--v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
xlabel('J = V/nD','FontWeight','bold')
ylabel('\eta','FontWeight','bold')
lg = plot(NaN,NaN,'k-d',NaN,NaN,'k-^',NaN,NaN,'k-v',NaN,NaN,NaN,NaN,'--');
for i = 1:length(lg)
    lg(i).LineWidth = 1.2;
end
lg(4).Color = [0 0.4470 0.7410];
lg(5).Color = [0.8500 0.3250 0.0980];
legend(lg,'Crigler Dual','Davidson Dual','Liebeck Dual','Davidson Analysis','Liebeck Analysis','Location','best','FontWeight','bold','Color','none')
axis([0.5 2.6 0.4 1])
grid minor; grid on

subplot(3,3,7)
title('g)','FontWeight','bold')
hold on
plot(ddcd.J,ddcd.eta,'-d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(dddd.J,dddd.eta,'-^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ddld.J,ddld.eta,'-v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ldcd.J,ldcd.eta,'--d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(lddd.J,lddd.eta,'--^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(ldld.J,ldld.eta,'--v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])

plot(lscs.J,lscs.eta,':o','LineWidth',1.5,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(lsls.J,lsls.eta,':s','LineWidth',1.5,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(cscs.J,cscs.eta,'-.o','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
plot(csls.J,csls.eta,'-.s','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
xlabel('J = V/nD','FontWeight','bold')
ylabel('\eta','FontWeight','bold')
axis([0.5 2.6 0.4 1])
grid minor; grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
subplot(3,3,2)
title('b)','FontWeight','bold')
hold on
plot(lscs.J,lscs.Ct,':o','LineWidth',1.5,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(lsls.J,lsls.Ct,':s','LineWidth',1.5,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(cscs.J,cscs.Ct,'-.o','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
plot(csls.J,csls.Ct,'-.s','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
xlabel('J = V/nD','FontWeight','bold')
ylabel('C_T','FontWeight','bold')
axis([0.5 2.6 0 0.45])
lg = plot(NaN,NaN,'k-o',NaN,NaN,'k-s',NaN,NaN,':',NaN,NaN,'-.');
for i = 1:length(lg)
    lg(i).LineWidth = 1.2;
end
lg(3).Color = [0.9290 0.6940 0.1250];
lg(4).Color =[0.4940 0.1840 0.5560];
legend(lg,'Crigler Single','Liebeck Single','Liebeck Analysis','Davidson Mod. Analysis','Location','best','FontWeight','bold','Color','none')
grid minor; grid on

subplot(3,3,5)
title('e)','FontWeight','bold')
hold on
plot(ddcd.J,ddcd.Ct,'-d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(dddd.J,dddd.Ct,'-^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ddld.J,ddld.Ct,'-v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ldcd.J,ldcd.Ct,'--d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(lddd.J,lddd.Ct,'--^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(ldld.J,ldld.Ct,'--v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
xlabel('J = V/nD','FontWeight','bold')
ylabel('C_T','FontWeight','bold')
axis([0.5 2.6 0 0.45])
lg = plot(NaN,NaN,'k-d',NaN,NaN,'k-^',NaN,NaN,'k-v',NaN,NaN,NaN,NaN,'--');
for i = 1:length(lg)
    lg(i).LineWidth = 1.2;
end
lg(4).Color = [0 0.4470 0.7410];
lg(5).Color = [0.8500 0.3250 0.0980];
legend(lg,'Crigler Dual','Davidson Dual','Liebeck Dual','Davidson Analysis','Liebeck Analysis','Location','best','FontWeight','bold','Color','none')
grid minor; grid on


subplot(3,3,8)
title('h)','FontWeight','bold')
hold on
plot(ddcd.J,ddcd.Ct,'-d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(dddd.J,dddd.Ct,'-^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ddld.J,ddld.Ct,'-v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ldcd.J,ldcd.Ct,'--d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(lddd.J,lddd.Ct,'--^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(ldld.J,ldld.Ct,'--v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])

plot(lscs.J,lscs.Ct,':o','LineWidth',1.5,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(lsls.J,lsls.Ct,':s','LineWidth',1.5,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(cscs.J,cscs.Ct,'-.o','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
plot(csls.J,csls.Ct,'-.s','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
xlabel('J = V/nD','FontWeight','bold')
ylabel('C_T','FontWeight','bold')
axis([0.5 2.6 0 0.45])
grid minor; grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,3,3)
title('c)','FontWeight','bold')
hold on
plot(lscs.J,lscs.Cp,':o','LineWidth',1.5,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(lsls.J,lsls.Cp,':s','LineWidth',1.5,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(cscs.J,cscs.Cp,'-.o','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
plot(csls.J,csls.Cp,'-.s','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
xlabel('J = V/nD','FontWeight','bold')
ylabel('C_P','FontWeight','bold')
axis([0.5 2.6 0 0.7])
lg = plot(NaN,NaN,'k-o',NaN,NaN,'k-s',NaN,NaN,':',NaN,NaN,'-.');
for i = 1:length(lg)
    lg(i).LineWidth = 1.2;
    
end
lg(3).Color = [0.9290 0.6940 0.1250];
lg(4).Color =[0.4940 0.1840 0.5560];
legend(lg,'Crigler Single','Liebeck Single','Liebeck Analysis','Davidson Mod. Analysis','Location','best','FontWeight','bold','Color','none')
grid minor; grid on

subplot(3,3,6)
title('f)','FontWeight','bold')
hold on
plot(ddcd.J,ddcd.Cp,'-d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(dddd.J,dddd.Cp,'-^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ddld.J,ddld.Cp,'-v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ldcd.J,ldcd.Cp,'--d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(lddd.J,lddd.Cp,'--^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(ldld.J,ldld.Cp,'--v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
xlabel('J = V/nD','FontWeight','bold')
ylabel('C_P','FontWeight','bold')
axis([0.5 2.6 0 0.7])
lg = plot(NaN,NaN,'k-d',NaN,NaN,'k-^',NaN,NaN,'k-v',NaN,NaN,NaN,NaN,'--');
for i = 1:length(lg)
    lg(i).LineWidth = 1.2;
end
lg(4).Color = [0 0.4470 0.7410];
lg(5).Color = [0.8500 0.3250 0.0980];
legend(lg,'Crigler Dual','Davidson Dual','Liebeck Dual','Davidson Analysis','Liebeck Analysis','Location','best','FontWeight','bold','Color','none')
grid minor; grid on


subplot(3,3,9)
title('i)','FontWeight','bold')
hold on
plot(ddcd.J,ddcd.Cp,'-d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(dddd.J,dddd.Cp,'-^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ddld.J,ddld.Cp,'-v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0 0.4470 0.7410])
plot(ldcd.J,ldcd.Cp,'--d','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(lddd.J,lddd.Cp,'--^','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])
plot(ldld.J,ldld.Cp,'--v','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.8500 0.3250 0.0980])

plot(lscs.J,lscs.Cp,':o','LineWidth',1.5,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(lsls.J,lsls.Cp,':s','LineWidth',1.5,'MarkerIndices',1:3:27,'Color',[0.9290 0.6940 0.1250])
plot(cscs.J,cscs.Cp,'-.o','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
plot(csls.J,csls.Cp,'-.s','LineWidth',1,'MarkerIndices',1:3:27,'Color',[0.4940 0.1840 0.5560])
xlabel('J = V/nD','FontWeight','bold')
ylabel('C_P','FontWeight','bold')
axis([0.5 2.6 0 0.7])
% legend('Liebeck Single','Crigler Single')
grid minor; grid on

%% Polares

load('polares3602.mat')
load('polares3602.mat')

figure
subplot(1,2,1)
hold on
plot(p.alpha,p.CL,'LineWidth',2)
plot(pol.alpha,pol.CL,'--')
xlabel('\alpha [°]','FontWeight','bold')
ylabel('C_L','FontWeight','bold')
legend('XFoil','Extrapolada','Location','best','FontWeight','bold','Color','none')
axis([-10 90 -0.2 1.6])
grid on; grid minor;

subplot(1,2,2)
hold on
plot(p.alpha,p.CD,'LineWidth',2)
plot(pol.alpha,pol.CD,'--')
xlabel('\alpha [°]','FontWeight','bold')
ylabel('C_D','FontWeight','bold')
legend('XFoil','Extrapolada','Location','best','FontWeight','bold','Color','none')
axis([-10 90 0 2])
grid on; grid minor;
