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
%load('Workspaces/Calibration','s_bar','s_max','verd','S_bar','sig','gamma','S');
load('WS_Krogh_Jrensen_Calibration_PD_Claim_1','S')
load('WS_Krogh_Jrensen_Calibration_PD_Claim_1','s_bar','s_max','verd','S_bar','sig','gamma','S','astsim');
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
Rec_s_bar = fzero(@(x) (integral(@q_s,-Inf,x) - rec_emp_percentage), s_bar-0.1);
%Rec_s_bar = -2.22;

[heights location] = hist(astsim, 500);
width = location(2) - location(1);
heights = heights / (size(astsim, 1) * width);

warning('off','all'); % fplot doest like the integral functions
figure;
bar(location, heights,'hist')
hold on
fplot(@q_s, [min(log(S)+3) s_max+0.15],'red');title('Stationary Distribution of s');
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