addpath('Functions');
addpath('Data');
addpath('Workspaces');
addpath('Calibration');
addpath('Figures');
addpath('Tables');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loads simulated data and performs regressions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
%%
PD_Claim_Regressions = 0; % 0 = PC
                          % 1 = PD

momPC = readtable('PC_Claim_Sim_mom.txt');
momPD = readtable('PD_Claim_Sim_mom.txt');
datPC = readtable('PC_Claim_Sim_dat.txt');
datPD = readtable('PD_Claim_Sim_dat.txt');
figure;
subplot(2,1,1)
plot(table2array(datPC(:,4)));title('$P/C$','Interpreter','latex')
subplot(2,1,2)
plot(table2array(datPD(:,4)));title('$P/D$','Interpreter','latex')
%%
if PD_Claim_Regressions == 0
    Moments = momPC;
    Data = datPC;
else
    Moments = momPD;
    Data = datPD;
end

aretssim =  table2array(Data(:,4));
astsim = table2array(Data(:,1));
if PD_Claim_Regressions == 0
    load('PC_Claim_workspace','s_bar','s_max',...
        'verd','S_bar','sig','gamma','S','astsim_pf','alnrtsim_pf');
else
    load('PD_Claim_workspace','s_bar','s_max',...
        'verd','S_bar','sig','gamma','S','astsim_pf','alnrtsim_pf');
end

%% Matching the empirical density
NBER_REC = importdata('USREC.csv');
% Define period yyyy-mm-dd
from = '1950-01-01';
to   = '2018-12-01';

% find indexes
idx_from = find(NBER_REC.textdata(:,1)==string(from)) - 1;
idx_to   = find(NBER_REC.textdata(:,1)==string(to)) - 1;

% Calculate percentage of the time the economy is in recession
rec_emp_percentage = sum(NBER_REC.data(idx_from:idx_to,1)) / length(NBER_REC.data(idx_from:idx_to,1));
rec_emp_percentagetotal = sum(NBER_REC.data(1:end,1)) / length(NBER_REC.data(1:end,1));
%%
s_bar_2 = log(0.02); % Recession specification by figure
Rec_s_bar = fzero(@(x) (integral(@q_s,-Inf,x) - rec_emp_percentage), s_bar-0.1);
Model_Rec = integral(@q_s,-Inf,s_bar);
Model_Rec_2 = integral(@q_s,-Inf,s_bar_2);
Match_Rec = integral(@q_s,-Inf,Rec_s_bar);
%Rec_s_bar = -2.22;
%%
load('PC_Claim_workspace','astsim_pf');astsim = astsim_pf;
[heights location] = hist(astsim, 75);
width = location(2) - location(1);
heights = heights / (size(astsim, 1) * width);
%%
warning('off','all'); % fplot doesnt like the integral functions
figure;
barplot = bar(location, heights,'hist');
barplot.FaceColor = [0, 0.4470, 0.7410];
hold on
fplot(@q_s, [min(log(S)-0.5) s_max+0.25],'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2.5);%title('Stationary Distribution of s');
hold on
xline(Rec_s_bar,'--','$\bar{s}_{rec}$','Interpreter','latex','FontSize',18);
hold on
xline(s_bar,'--','$\bar{s}$','Interpreter','latex','FontSize',18);
xline(s_bar_2,'--','$\bar{s}_{2,rec}$','Interpreter','latex','FontSize',18);
xlim([min(log(S))-0.5 -2]);
legend('Histogram','Theoretical Density','Location','northwest')
hold off
saveas(gcf,'../Figures/DistributionS_t','epsc')
%% Redefining recession periods in the simulation
% such that the frequency of recession in the simulation corresponds to the
% empirical frequency of recessions:
% Recession s_t < Rec_s_bar
rec_sim_ss = NaN(length(astsim), 1);
for i = 1:length(astsim)
    if astsim(i) < Rec_s_bar
        rec_sim_ss(i) = 1;
    else
        rec_sim_ss(i) = 0;
    end
end
rec_sim_ss_percentage = sum(rec_sim_ss(:)==1) / length(rec_sim_ss);
%%
load('PD_Claim_workspace','s_bar','s_max',...
        'verd','S_bar','sig','gamma','S','astsim','alnrtsim_pf','alnpctsim_pf','Erfinterp_pf');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% r_(t+h) = alpha + beta p/d_t + eps %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PD_regress = alnpctsim_pf;              % PD/PC
rfr  = Erfinterp_pf;                    % Risk free rate
rets = alnrtsim_pf - rfr;               % Excess Returns
h    =  1;                              % Forecast Horizon 0 = in-sample regression

y   = rets(1+h:end,1);
x   = [ones(length(rets(1:end-h,:)), 1),  ...            
    rec_sim_ss(1:end-h,:) .* PD_regress(1:end-h,1), ...  
    (1-rec_sim_ss(1:end-h,:)) .* PD_regress(1:end-h,1)]; 
regPDrec = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
    rec_sim_ss(1:end-h,:) .* PD_regress(1:end-h,1)]; 
regPDrec1 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
    (1-rec_sim_ss(1:end-h,:)) .* PD_regress(1:end-h,1)]; 
regPDexp1 = nwest(y,x,0);

% no business cycle
x   = [ones(length(rets(1:end-h,:)), 1),...  
    PD_regress(1:end-h,1)];
regPDnorec = nwest(y,x,0);

load('PC_Claim_workspace','alnrtsim_pf','alnpctsim_pf','Erfinterp_pf');
PC_regress = alnpctsim_pf;              % PD/PC
rfr  = Erfinterp_pf;                    % Risk free rate
rets = alnrtsim_pf - rfr;               % Excess Returns

y   = rets(1+h:end,1);

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
    rec_sim_ss(1:end-h,:) .* PC_regress(1:end-h,1), ...  
    (1-rec_sim_ss(1:end-h,:)) .* PC_regress(1:end-h,1)]; 
regPCrec = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
    rec_sim_ss(1:end-h,:) .* PC_regress(1:end-h,1)]; 
regPCrec1 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
    (1-rec_sim_ss(1:end-h,:)) .* PC_regress(1:end-h,1)]; 
regPCexp1 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...       
    PC_regress(1:end-h,1)];
regPCnorec = nwest(y,x,0);
regs1 = [regPCrec regPDrec regPCrec1 regPCexp1 regPDrec1 regPDexp1];
% regs = [regPDrec regPDnorec regPCrec regPCnorec];
RegressionTable2;
%% Regressions 2
for i = 1:length(astsim)
    if astsim(i) < log(0.02)
        rec_sim_02(i) = 1;
    else
        rec_sim_02(i) = 0;
    end 
end
rec_sim_02 = rec_sim_02';
%%
y   = rets(1+h:end,1);
x   = [ones(length(rets(1:end-h,:)), 1),  ...         
    rec_sim_02(1:end-h,:) .* PC_regress(1:end-h,1), ...
    (1-rec_sim_02(1:end-h,:)) .* PC_regress(1:end-h,1)];
regPCrec2 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...            
    rec_sim_02(1:end-h,:) .* PD_regress(1:end-h,1), ...  
    (1-rec_sim_02(1:end-h,:)) .* PD_regress(1:end-h,1)]; 
regPDrec2 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...        
    (1-rec_sim_02(1:end-h,:)) .* PC_regress(1:end-h,1)];
regPCexp3 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
     rec_sim_02(1:end-h,:) .* PC_regress(1:end-h,1)];
regPCrec3 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
     rec_sim_02(1:end-h,:) .* PD_regress(1:end-h,1)];
regPDrec3 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
     (1 - rec_sim_02(1:end-h,:)) .* PD_regress(1:end-h,1)];
regPDexp3 = nwest(y,x,0);

regs2 = [regPCrec2 regPDrec2 regPCrec3 regPCexp3 regPDrec3 regPDexp3];
RegressionTable1;
%% Split
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regime Switching model Observable states %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('PC_Claim_workspace','alnrtsim_pf','alnpctsim_pf','Erfinterp_pf');
PC_regress = alnpctsim_pf;              % PD/PC
rfr  = Erfinterp_pf;                    % Risk free rate
rets = alnrtsim_pf - rfr;               % Excess Returns

yrec  = rec_sim_ss(1+h:end,:) .* rets(1+h:end,1); % Recession
yexp  = (1-rec_sim_ss(1+h:end,:)) .* rets(1+h:end,1); 
xrec  = [ones(length(rets(1:end-h,:)), 1),  ...         
        rec_sim_ss(1:end-h,:) .* PC_regress(1:end-h,1)];
xexp  = [ones(length(rets(1:end-h,:)), 1),...
        (1-rec_sim_ss(1:end-h,:)) .* PC_regress(1:end-h,1)];
RegRec = nwest(yrec,xrec,0);
RegExp = nwest(yexp,xexp,0);
%% Split
load('PD_Claim_workspace','alnrtsim_pf','alnpctsim_pf','Erfinterp_pf');
PC_regress = alnpctsim_pf;              % PD/PC
rfr  = Erfinterp_pf;                    % Risk free rate
rets = alnrtsim_pf - rfr;               % Excess Returns

yrec  = rec_sim_ss(1+h:end,:) .* rets(1+h:end,1); % Recession
yexp  = (1-rec_sim_ss(1+h:end,:)) .* rets(1+h:end,1); 
xrec  = [ones(length(rets(1:end-h,:)), 1),...         
        rec_sim_ss(1:end-h,:) .* PC_regress(1:end-h,1)];
xexp  = [ones(length(rets(1:end-h,:)), 1),...
        (1-rec_sim_ss(1:end-h,:)) .* PC_regress(1:end-h,1)];
RegRec_PD = nwest(yrec,xrec,0);
RegExp_PD = nwest(yexp,xexp,0);
RSregs = [RegRec, RegRec_PD, RegExp, RegExp_PD];
RSRegressionTable;
%%
load('PD_Claim_workspace','elnrtsim'); ExpRetsPD = elnrtsim;
load('PC_Claim_workspace','elnrtsim'); ExpRetsPC = elnrtsim;
%%
figure;
subplot(2,1,1)
plot(ExpRetsPC);title({'$P/C$', ['$E( E_t  (r_{t+1}) )$ =',num2str(mean(ExpRetsPC),6)]},'Interpreter','latex');
xlim([-500 100000])
ylabel('$ E_t  (r_{t+1})$','FontSize',14,'interpreter','latex');
subplot(2,1,2)
plot(ExpRetsPD);title({'$P/D$', ['$E( E_t  (r_{t+1}) )$ =',num2str(mean(ExpRetsPD),6)]},'Interpreter','latex');
ylabel('$ E_t  (r_{t+1})$','FontSize',14,'interpreter','latex');
xlim([-500 100000]);
saveas(gcf,'../Figures/Excess_Rets','epsc');
%% Persistence s_t
x = astsim_pf(1:end-1,:);
x1 = astsim_pf(2:end,:);
autocorrX = x\x1
%%
load('PD_Claim_workspace','alnpctsim_pf'); PDratio = alnpctsim_pf;
load('PC_Claim_workspace','alnpctsim_pf'); PCratio = alnpctsim_pf;
subplot(2,1,1)
plot(PCratio);ylabel('$p_t-c_t$','FontSize',14,'Interpreter','latex');
title({'$P/C$', ['$E(p_t-c_t)$ =',num2str(mean(PCratio),4)]},'Interpreter','latex');
subplot(2,1,2)
plot(PDratio);ylabel('$p_t-d_t$','FontSize',14,'Interpreter','latex');
title({'$P/D$', ['$E(p_t-d_t)$ =',num2str(mean(PDratio),4)]},'Interpreter','latex');
ylim([1.25 3.5]);
saveas(gcf,'../Figures/PCPD_chain','epsc');
%%
clear
load('PD_Claim_workspace','Edc_pf', 'Stdc_pf', 'Erfinterp_pf', 'Shprinterp_pf', 'ShpRinterp_pf', 'Eexrettinterp_pf', 'Stdexrettinterp_pf', 'Ep_d_pf', 'Stdp_d_pf');
PD_edc = Edc_pf;
PD_stdc = Stdc_pf;
PD_rfr = Erfinterp_pf;
PD_shrp = Shprinterp_pf;
PD_SHRP = ShpRinterp_pf;
PD_Eer = Eexrettinterp_pf;
PD_Stder = Stdexrettinterp_pf;
PD_Epd = Ep_d_pf;
PD_StdPD = Stdp_d_pf;
load('PC_Claim_workspace','Edc_pf', 'Stdc_pf', 'Erfinterp_pf', 'Shprinterp_pf', 'ShpRinterp_pf', 'Eexrettinterp_pf', 'Stdexrettinterp_pf', 'Ep_d_pf', 'Stdp_d_pf');
PC_edc = Edc_pf;
PC_stdc = Stdc_pf;
PC_rfr = Erfinterp_pf;
PC_shrp = Shprinterp_pf;
PC_SHRP = ShpRinterp_pf;
PC_Eer = Eexrettinterp_pf;
PC_Stder = Stdexrettinterp_pf;
PC_Epd = Ep_d_pf;
PC_StdPD = Stdp_d_pf;
load('CC_PC_Claim_workspace','Edc_pf', 'Stdc_pf', 'Erfinterp_pf', 'Shprinterp_pf', 'ShpRinterp_pf', 'Eexrettinterp_pf', 'Stdexrettinterp_pf', 'Ep_d_pf', 'Stdp_d_pf');
CC_PC_edc = Edc_pf;
CC_PC_stdc = Stdc_pf;
CC_PC_rfr = Erfinterp_pf;
CC_PC_shrp = Shprinterp_pf;
CC_PC_SHRP = ShpRinterp_pf;
CC_PC_Eer = Eexrettinterp_pf;
CC_PC_Stder = Stdexrettinterp_pf;
CC_PC_Epd = Ep_d_pf;
CC_PC_StdPD = Stdp_d_pf;
load('CC_PD_Claim_workspace','Edc_pf', 'Stdc_pf', 'Erfinterp_pf', 'Shprinterp_pf', 'ShpRinterp_pf', 'Eexrettinterp_pf', 'Stdexrettinterp_pf', 'Ep_d_pf', 'Stdp_d_pf');
CC_PD_edc = Edc_pf;
CC_PD_stdc = Stdc_pf;
CC_PD_rfr = Erfinterp_pf;
CC_PD_shrp = Shprinterp_pf;
CC_PD_SHRP = ShpRinterp_pf;
CC_PD_Eer = Eexrettinterp_pf;
CC_PD_Stder = Stdexrettinterp_pf;
CC_PD_Epd = Ep_d_pf;
CC_PD_StdPD = Stdp_d_pf;
Simulatedmom;
