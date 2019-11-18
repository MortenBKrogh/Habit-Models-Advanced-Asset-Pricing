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
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
%%
PD_Claim_Regressions = 1; % 0 = PC
                          % 1 = PD

momPC = readtable('PC_Claim_Sim_mom.txt');
momPD = readtable('PD_Claim_Sim_mom.txt');
datPC = readtable('PC_Claim_Sim_dat.txt');
datPD = readtable('PD_Claim_Sim_dat.txt');
figure;
subplot(2,1,1)
plot(table2array(datPC(:,4)));title('$P/C$')
subplot(2,1,2)
plot(table2array(datPD(:,4)));title('$P/D$')
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
        'verd','S_bar','sig','gamma','S','astsim','alnrtsim_pf');
else
    load('PD_Claim_workspace','s_bar','s_max',...
        'verd','S_bar','sig','gamma','S','astsim','alnrtsim_pf');
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
gamma = 2
Rec_s_bar = fzero(@(x) (integral(@q_s,-Inf,x) - rec_emp_percentage), s_bar-0.1);
Model_Rec = integral(@q_s,-Inf,s_bar);
Match_Rec = integral(@q_s,-Inf,Rec_s_bar);
%Rec_s_bar = -2.22;
%%
[heights location] = hist(astsim, 65);
width = location(2) - location(1);
heights = heights / (size(astsim, 1) * width);
%%
warning('off','all'); % fplot doesnt like the integral functions
figure;
barplot = bar(location, heights,'hist');
barplot.FaceColor = [0, 0.4470, 0.7410];
hold on
fplot(@q_s, [min(log(S)+1.5) s_max+0.15],'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2.5);%title('Stationary Distribution of s');
hold on
xline(Rec_s_bar,'--','$\bar{s}_{rec}$','Interpreter','latex','FontSize',18);
hold on
xline(s_bar,'--','$\bar{s}$','Interpreter','latex','FontSize',18);
xlim([min(log(S)+1.5) -2]);
legend('Histogram','Theoretical Density','Location','northwest')
hold off
%%saveas(gcf,'../Figures/DistributionS_t','epsc')
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
x   = [ones(length(rets(1:end-h,:)), 1),  ...            % vector of Ones
    rec_sim_ss(1:end-h,:) .* PD_regress(1:end-h,1), ...  %    I_rec_t *PD_t
    (1-rec_sim_ss(1:end-h,:)) .* PD_regress(1:end-h,1)]; % (1-I_rec_t)*PD_t
regPDrec = nwest(y,x,0);

%%
y   = rets(1+h:end,1);
x   = [ones(length(rets(1:end-h,:)), 1),...         % vector of Ones
    PD_regress(1:end-h,1)];
regPDnorec = nwest(y,x,0);
%%
load('PC_Claim_workspace','alnrtsim_pf','alnpctsim_pf','Erfinterp_pf');
PC_regress = alnpctsim_pf;              % PD/PC
rfr  = Erfinterp_pf;                    % Risk free rate
rets = alnrtsim_pf - rfr;               % Excess Returns

y   = rets(1+h:end,1);
x   = [ones(length(rets(1:end-h,:)), 1),  ...         % vector of Ones
    rec_sim_ss(1:end-h,:) .* PC_regress(1:end-h,1), ...  %    I_rec_t *PD_t
    (1-rec_sim_ss(1:end-h,:)) .* PC_regress(1:end-h,1)]; % (1-I_rec_t)*PD_t
regPCrec = nwest(y,x,0);
%%
y   = rets(1+h:end,1);
x   = [ones(length(rets(1:end-h,:)), 1),  ...         % vector of Ones
    PC_regress(1:end-h,1)];
regPCnorec = nwest(y,x,0);
regs = [regPDrec regPDnorec regPCrec regPCnorec];
RegressionTable;
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
xrec  = [ones(length(rets(1:end-h,:)), 1),  ...         
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
plot(ExpRetsPC);title({'$P/C$', ['mean =',num2str(mean(ExpRetsPC),6)]});
xlim([-500 100000])
ylabel('Excess Returns');
subplot(2,1,2)
plot(ExpRetsPD);title({'$P/D$', ['mean =',num2str(mean(ExpRetsPD),6)]});
ylabel('Excess Returns');
xlim([-500 100000]);
% saveas(gcf,'../Figures/Excess_Rets','epsc');
%% Persistence s_t
x = astsim_pf(1:end-1,:);
x1 = astsim_pf(2:end,:);
autocorrX = x\x1