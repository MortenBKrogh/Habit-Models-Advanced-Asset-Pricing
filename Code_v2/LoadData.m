%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUNNING THIS FILE PRODUCES MOST THE TABLES AND FIGURES IN THE ARTICLE 
% JENSEN & KROGH (2019).
% The program relies on running the MAIN.m file according to the
% description in that file. This is because this program is working upon
% the 4 workspaces saved when running the MAIN.m program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
opts.Colors     = get(groot,'defaultAxesColorOrder');

%%
Save_Figures = 0;         % 0 = dont save
                          % 1 = save    
                          
load('PC_Claim_workspace','s_bar','s_max',...
         'verd','S_bar','sig','gamma','S','astsim_pf','alnrtsim_pf','stsim');

%% Matching the empirical receession probability with the models density of s_t
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

labeledMatrix = bwlabel(NBER_REC.data(idx_from:idx_to,1));
measurements = regionprops(labeledMatrix, 'Area');
RecessionLengths = [measurements.Area];
mean(RecessionLengths)/12;
%%
s_bar_2 = log(0.02); % Recession specification by figure
Rec_s_bar = fzero(@(x) (integral(@q_s,-Inf,x) - rec_emp_percentage), s_bar-0.1);
Model_Rec = integral(@q_s,-Inf,s_bar);
Model_Rec_2 = integral(@q_s,-Inf,s_bar_2);
Match_Rec = integral(@q_s,-Inf,Rec_s_bar);
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
if Save_Figures
saveas(gcf,'../Figures/DistributionS_t','epsc')
end
%% Risk aversion plot
RA = gamma./exp(stsim);
mRA = mean(RA);
figure;
plot(RA);ylabel('$\gamma/S_t,\qquad$ $\gamma=2$','Interpreter','latex');
xlabel('$t$','interpreter','latex');
yline(mean(RA),'--','LineWidth',2);
legend('$\gamma/S_t$: Risk Aversion',['$\mathbf{E}\gamma/S_t$=',num2str(mRA)],'interpreter','latex')
xlim([0 100000]);
if Save_Figures
saveas(gcf,'../Figures/RA','epsc')
end
max(RA)
%% Redefining recession periods in the simulation
% such that the frequency of recession in the simulation corresponds to the
% empirical frequency of recessions:
% Recession s_t < Rec_s_bar
load('PC_Claim_workspace','stsim');astsim = stsim;
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
labeledMatrix = bwlabel(rec_sim_ss);
measurements = regionprops(labeledMatrix, 'Area');
RecessionLengths = [measurements.Area];
mean(RecessionLengths)
%%
    load('PD_Claim_workspace','s_bar','s_max',...
        'verd','S_bar','sig','gamma','S','stsim','lnrtsim','lnpctsim','Erfinterp_pf');
    Erfinterp_pf = Erfinterp_pf./12;
    PD_regress = lnpctsim(2:end,1);              % PD
    lnrtsim_PD = lnrtsim;
    load('PC_Claim_workspace','lnrtsim','lnpctsim')
    lnrtsim_PC = lnrtsim;
    PC_regress = lnpctsim(2:end,1);              % PC
    rec_sim_ss = rec_sim_ss(2:end,1);
%%
rfr  = Erfinterp_pf;                       % Risk free rate
Erets_PC = lnrtsim_PC - rfr;               % Excess Returns PC
Erets_PD = lnrtsim_PD - rfr;
h    =  1;                                   % Forecast Horizon 0 = in-sample regression
yPC  = Erets_PC(1+h:end,1);                  % Regressand PC
yPD  = Erets_PD(1+h:end,1);                  % Regressand PD
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    No business cycle regressions   %%%
%%% r_(t+h) = alpha + beta p/d_t + eps %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x   = [ones(length(PD_regress(1:end-h,1)), 1),...  
    PD_regress(1:end-h,1)];
regPDnorec = nwest(yPD,x,0);

x   = [ones(length(PC_regress(1:end-h,1)), 1),  ...       
    PC_regress(1:end-h,1)];
regPCnorec = nwest(yPC,x,0);

regsNB = [regPCnorec regPDnorec];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Business cycle regressions                   %%%
%%% r_(t+h) = alpha + beta_1 p/d_t*I_rec + beta_2(1-I_rec)p/d_t + eps %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Full Business cycle                         %
x   = [ones(length(PD_regress(1:end-h,1)), 1),  ...            
    rec_sim_ss(1:end-h,:) .* PD_regress(1:end-h,1), ...  
    (1-rec_sim_ss(1:end-h,:)) .* PD_regress(1:end-h,1)]; 
regPDrec = nwest(yPD,x,0); %% Full BC <- PD

x   = [ones(length(PC_regress(1:end-h,1)), 1),  ...         
    rec_sim_ss(1:end-h,:) .* PC_regress(1:end-h,1), ...  
    (1-rec_sim_ss(1:end-h,:)) .* PC_regress(1:end-h,1)]; 
regPCrec = nwest(yPC,x,0); %% Full BC <- PC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Split Business cycle                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lower_Sbar = 0; %% Table purposes only do not change

retsHRecPD = yPD .* rec_sim_ss(2:end);     %% Excess Returns Recession
retsHExpPD = yPD .* (1-rec_sim_ss(2:end)); %% Exceess Returns Expansions
retsHRecPC = yPC .* rec_sim_ss(2:end);     %% Excess Returns Recession
retsHExpPC = yPC .* (1-rec_sim_ss(2:end)); %% Exceess Returns Expansions

PDRegHRec = rec_sim_ss(1:end-h,:) .* PD_regress(1:end-h,1);
PDRegHExp =  (1 - rec_sim_ss(1:end-h,:)) .* PD_regress(1:end-h,1);
PCRegHRec = rec_sim_ss(1:end-h,:) .* PD_regress(1:end-h,1);
PCRegHExp =  (1 - rec_sim_ss(1:end-h,:)) .* PC_regress(1:end-h,1);

a = [retsHRecPD, PDRegHRec];
a = a(all(a,2),:);
ExcRetsRec = a(:,1);                   %% <- Excess Returns Recessions only
ExRetsRecPDFC = ExcRetsRec;
PDrecHR = [ones(size(a,1), 1) a(:,2)]; %% <- PD recession
regPDrec1 = nwest(ExcRetsRec,PDrecHR,0); 

a = [retsHExpPD, PDRegHExp];
a = a(all(a,2),:);
ExRetsExp = a(:,1); %% <- Excess Returns Expansions only
ExRetsExpPDFC = ExRetsExp;
PDexpHR   = [ones(size(a,1),1) a(:,2)]; %% <- PD Expansion
regPDexp1 = nwest(ExRetsExp,PDexpHR,0);

a = [retsHExpPC, PCRegHExp];
a = a(all(a,2),:);
ExRetsExp = a(:,1);
ExRetsExpPCFC = ExRetsExp;
PCExpHR   = [ones(size(a,1),1) a(:,2)];
regPCexp1 = nwest(ExRetsExp,PCExpHR,0);

a = [retsHRecPC, PCRegHRec];
a = a(all(a,2),:);
ExRecRets = a(:,1);
ExRetsRecPCFC = ExRecRets;
PCrecHR   = [ones(size(a,1),1) a(:,2)];
regPCrec1 = nwest(ExRecRets,PCrecHR,0);

regs1 = [regPCrec regPDrec regPCrec1 regPCexp1 regPDrec1 regPDexp1];
if Save_Figures
RegressionTable2;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Split Business cycle lower s bar                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('PD_Claim_workspace','s_bar','s_max',...
        'verd','S_bar','sig','gamma','S','stsim','lnrtsim','lnpctsim','Erfinterp_pf');
    Erfinterp_pf = Erfinterp_pf./12;
    PD_regress   = lnpctsim(2:end,1);             % PD
    lnrtsimPD    = lnrtsim;
    load('PC_Claim_workspace','lnrtsim','lnpctsim')
    alnrtsim_pf = lnrtsim;
    PC_regress  = lnpctsim(2:end,1);              % PC
    lnrtsimPC   = lnrtsim; 
    h=1;
rfr  = Erfinterp_pf;                    % Risk free rate
retsPC = lnrtsimPC - rfr;               % Excess Returns
retsPD = lnrtsimPD - rfr;
h    =  1;                              % Forecast Horizon 0 = in-sample regression
yPC  =  retsPC(1+h:end,1);              % Regressand 
yPD  =  retsPD(1+h:end,1);
rec_sim_02 = zeros(size(stsim,1),1);
for i = 1:size(stsim,1)
    if stsim(i) < log(0.02)
        rec_sim_02(i) = 1;
    else
        rec_sim_02(i) = 0;
    end 
end
rec_sim_02 = rec_sim_02(2:end);
%%
lower_Sbar = 1;
s_bar_2 = log(0.02);
x   = [ones(length(yPD), 1),  ...            
    rec_sim_02(1:end-h,:) .* PD_regress(1:end-h,1), ...  
    (1-rec_sim_02(1:end-h,:)) .* PD_regress(1:end-h,1)]; 
regPDrec = nwest(yPD,x,0); %% Full BC <- PD

x   = [ones(length(yPC), 1),  ...         
    rec_sim_02(1:end-h,:) .* PC_regress(1:end-h,1), ...  
    (1-rec_sim_02(1:end-h,:)) .* PC_regress(1:end-h,1)]; 
regPCrec = nwest(yPC,x,0); %% Full BC <- PC


retsHRecPC = retsPC(1+h:end) .* rec_sim_02(1+h:end);     %% Excess Returns Recession
retsHExpPC = retsPC(1+h:end) .* (1-rec_sim_02(1+h:end)); %% Excess Returns Expansions
retsHRecPD = retsPD(1+h:end) .* rec_sim_02(1+h:end);     %% Excess Returns Recession
retsHExpPD = retsPD(1+h:end) .* (1-rec_sim_02(1+h:end)); %% Excess Returns Expansions

PDRegHRec = rec_sim_02(1:end-h,:) .* PD_regress(1:end-h,1);
PDRegHExp =  (1 - rec_sim_02(1:end-h,:)) .* PD_regress(1:end-h,1);
PCRegHRec = rec_sim_02(1:end-h,:) .* PD_regress(1:end-h,1);
PCRegHExp =  (1 - rec_sim_02(1:end-h,:)) .* PC_regress(1:end-h,1);

a = [retsHRecPD, PDRegHRec];
a = a(all(a,2),:);
ExcRetsRec = a(:,1);                   %% <- Excess Returns Recessions only
PDrecHR1 = [ones(size(a,1), 1) a(:,2)]; %% <- PD recession
regPDrec1 = nwest(ExcRetsRec,PDrecHR1,0); 

a = [retsHExpPD, PDRegHExp];
a = a(all(a,2),:);
ExRetsExp = a(:,1);                     %% <- Excess Returns Expansions only
PDexpHR1   = [ones(size(a,1),1) a(:,2)]; %% <- PD Expansion
regPDexp1 = nwest(ExRetsExp,PDexpHR1,0);

a = [retsHExpPC, PCRegHExp];
a = a(all(a,2),:);
ExRetsExp = a(:,1);
PCExpHR1   = [ones(size(a,1),1) a(:,2)];
regPCexp1 = nwest(ExRetsExp,PCExpHR1,0);

a = [retsHRecPC, PCRegHRec];
a = a(all(a,2),:);
ExRecRets = a(:,1);
PCrecHR1   = [ones(size(a,1),1) a(:,2)];
regPCrec1 = nwest(ExRecRets,PCrecHR1,0);

regs2 = [regPCrec regPDrec regPCrec1 regPCexp1 regPDrec1 regPDexp1];
if Save_Figures
RegressionTable2;
end
%%
load('PD_Claim_workspace','elnrtsim','tsc'); ExpRetsPD = elnrtsim;
load('PC_Claim_workspace','elnrtsim'); ExpRetsPC = elnrtsim;
figure;
subplot(2,1,1)
plot(ExpRetsPC);title({'$P/C$', ['$E( E_t  (r_{t+1}) )$ =',num2str(mean(ExpRetsPC),6)]},'Interpreter','latex');
%xlim([-50 8333])
ylabel('$ E_t  (r_{t+1})$','FontSize',14,'interpreter','latex');
subplot(2,1,2)
plot(ExpRetsPD);title({'$P/D$', ['$E( E_t  (r_{t+1}) )$ =',num2str(mean(ExpRetsPD),6)]},'Interpreter','latex');
ylabel('$ E_t  (r_{t+1})$','FontSize',14,'interpreter','latex');
%xlim([-50 8333]);
if Save_Figures
saveas(gcf,'../Figures/Excess_Rets','epsc');
end
%% Persistence s_t
x = astsim_pf(1:end-1,:);
x1 = astsim_pf(2:end,:);
autocorrX = x\x1;
%%
load('PD_Claim_workspace','lnpctsim'); PDratio = lnpctsim;
load('PC_Claim_workspace','lnpctsim'); PCratio = lnpctsim;
name = '../Figures/PCPDMonthly_chain';
subplot(2,1,1)
plot(PCratio);ylabel('$p_t-c_t$','FontSize',14,'Interpreter','latex');
title({'$P/C$', ['$E(p_t-c_t)$ =',num2str(mean(PCratio),4)]},'Interpreter','latex');
xlim([1 100000]);
subplot(2,1,2)
plot(PDratio);ylabel('$p_t-d_t$','FontSize',14,'Interpreter','latex');
xlim([1 100000]);
title({'$P/D$', ['$E(p_t-d_t)$ =',num2str(mean(PDratio),4)]},'Interpreter','latex');
if Save_Figures
saveas(gcf,name,'epsc');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Long run regressions based on simulated data     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('PC_Claim_workspace','Erfinterp_pf','lnrtsim','lnpctsim');
rfr  = Erfinterp_pf;                 
retsPC = lnrtsim - (rfr/4); 
PCrat = lnpctsim(2:end);
load('PD_Claim_workspace','lnrtsim','lnpctsim')
PDrat = lnpctsim(2:end);
retsPD = lnrtsim - (rfr/4); 
j = 1;
T = length(PCrat);
Ta = T;
h= [1 2 3 5 7 10] * 12; 
while j <= size(h,2)
k = h(1,j);
xC = [ones(T-k+1,1), PCrat(1:T-k+1)];   
yC  = retsPC(1:T-k+1);%-rfr;
xD = [ones(T-k+1,1) PDrat(1:T-k+1)];
yD = retsPD(1:T-k+1);
    i = 2;
while i <= k
 yC = yC + retsPC(i:T-k+i);
 yD = yD + retsPD(i:T-k+i);
 i = i+1;
end   
b = xC\yC;
bmat(:,j) = b;
R2(:,j) = (std(xC * b)/ std(yC) )^2;
bd = xD\yD;
bdmat(:,j) = bd;
R2d(:,j) = (std(xD * bd)/ std(yD) )^2;
j = j+1;
end

tab = [bmat(2,:)', R2', bdmat(2,:)',  R2d']
names = split(char(num2str(h/12,1)));
varnames = split(['$\beta_{pc}$ ', '$R^2_{pc}$  ','$\beta_{pd}$ ','$R^2_{od}$'])'
tab =  array2table(tab,'Rownames',names,'VariableNames',varnames)
%% Moments Table
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
if Save_Figures
Simulatedmom;
end
%%
load('PC_Claim_workspace','stsim')
rec_sim_ss = NaN(length(stsim), 1);
for i = 1:length(stsim)
    if stsim(i) < Rec_s_bar
        rec_sim_ss(i) = 1;
    else
        rec_sim_ss(i) = 0;
    end
end
rec_sim_ss_percentage = sum(rec_sim_ss(:)==1) / length(rec_sim_ss);
%%
labeledMatrix = bwlabel(rec_sim_ss);
measurements = regionprops(labeledMatrix, 'Area');
RecessionLengths = [measurements.Area];
mean(RecessionLengths)/12;
%% FORECAST 1-period Expanding-Window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       Forecast Measures                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SampleSize = 0.05; % startng point 0 = full data [0:1[
WindowSize = 120;  % 240 = 20 years
%%%                          Recessions                                 %%%
init = SampleSize * size(ExRetsRecPDFC,1)+1;
MaxFC  = size(ExRetsRecPDFC,1)-WindowSize-1;
% LHS, Excess returns at time t+1
% RHS, ratios t
j = 1;
for i = init:MaxFC
    b1 = PDrecHR(i:i+WindowSize-1,:)\ExRetsRecPDFC(i:i+WindowSize-1,1);
    fit = b1(1) * PDrecHR(i+WindowSize,1) + b1(2) * PDrecHR(i+WindowSize,2);
    PDRecError(j,1) = ExRetsRecPDFC(i+WindowSize,:) - fit;
    PDRecMeanError(j,1) = ExRetsRecPDFC(i+WindowSize,:) - mean(ExRetsRecPDFC(i:i+WindowSize-1,:));
    
    b1 = PCrecHR(i:i+WindowSize-1,:)\ExRetsRecPCFC(i:i+WindowSize-1,1);
    fit = b1(1) * PCrecHR(i+WindowSize,1)+b1(2) * PCrecHR(i+WindowSize,2);
    PCRecError(j,1) = ExRetsRecPCFC(i+WindowSize,:) - fit;
    PCRecMeanError(j,1) = ExRetsRecPCFC(i+WindowSize,:) - mean(ExRetsRecPCFC(i:i+WindowSize-1,:));
    
    j=j+1;
end
%% Expansions
Init = SampleSize * size(ExRetsExpPDFC,1)+1;
MaxFC  = size(ExRetsExpPDFC,1)-WindowSize-1;
j = 1;
for i = init:MaxFC    
    b1 = PDexpHR(i:i+WindowSize-1,:)\ExRetsExpPDFC(i:i+WindowSize-1,1);
    fit = b1(1) * PDexpHR(i+WindowSize,1) + b1(2) * PDexpHR(i+WindowSize,2);
    PDExpError(j,1) = ExRetsExpPDFC(i+WindowSize,:) - fit;
    PDExpMeanError(j,1) = ExRetsExpPDFC(i+WindowSize,:) - mean(ExRetsExpPDFC(i:i+WindowSize-1,:));
    
    b1 = PCExpHR(i:i+WindowSize-1,:)\ExRetsExpPCFC(i:i+WindowSize-1,1);
    fit = b1(1) * PCExpHR(i+WindowSize,1) + b1(2) * PCExpHR(i+WindowSize,2);
    PCExpError(j,1) = ExRetsExpPCFC(i+WindowSize,:) - fit;
    PCExpMeanError(j,1) = ExRetsExpPCFC(i+WindowSize,:) - mean(ExRetsExpPCFC(i:i+WindowSize-1,:));
    
    j=j+1;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       R^2 OOS                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OoSR2ExpPD = 1-sum(PDExpError.^2)/sum(PDExpMeanError.^2);
OoSR2ExpPC = 1-sum(PCExpError.^2)/sum(PCExpMeanError.^2);
OoSR2recPD = 1-sum(PDRecError.^2)/sum(PDRecMeanError.^2);
OoSR2recPC = 1-sum(PCRecError.^2)/sum(PCRecMeanError.^2);
OoSR2 = [OoSR2recPC OoSR2recPD OoSR2ExpPC OoSR2ExpPD];