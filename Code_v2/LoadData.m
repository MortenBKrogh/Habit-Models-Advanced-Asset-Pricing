clear
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%%
momPC = readtable('PC_Claim_Sim_mom.txt');
momPD = readtable('PD_Claim_Sim_mom.txt');
datPC = readtable('PC_Claim_Sim_dat.txt');
datPD = readtable('PD_Claim_Sim_dat.txt');
%%
% table2latex(momPC, 'Tables/Moments_PC')
% table2latex(momPD, 'Tables/Moments_PD')
PDret = table2array(datPD(:,4));
PCret = table2array(datPC(:,4));

%%
astsim_pf = table2array(datPC(:,1));
%%
figure;
subplot(2,1,1)
plot(table2array(datPC(:,4)));title('$P/C$')
subplot(2,1,2)
plot(table2array(datPD(:,4)));title('$P/D$')
%%
NBER_REC = importdata('USREC.csv');

% Define period yyyy-mm-dd
from = '1950-01-01';
to   = '2018-12-01';

% find indexes
idx_from = find(NBER_REC.textdata(:,1)==string(from)) - 1;
idx_to   = find(NBER_REC.textdata(:,1)==string(to)) - 1;

% Calculate percentage of the time the economy is in recession
rec_emp_percentage = sum(NBER_REC.data(idx_from:idx_to,1)) / length(NBER_REC.data(idx_from:idx_to,1));
% define recession dummey as when s_t < s_bar
rec_sim_ss = NaN(length(astsim_pf), 1);
for i = 1:length(astsim_pf)
    if astsim_pf(i) < s_bar
        rec_sim_ss(i) = 1;
    else
        rec_sim_ss(i) = 0;
    end
end
rec_sim_ss_percentage = sum(rec_sim_ss(:)==1) / length(rec_sim_ss);
%% Matching the empirical density
Rec_s_bar = fzero(@(x) (integral(@q_s,-Inf,x) - rec_emp_percentage), s_bar-0.9);
%% Redefining recession periods in the simulation
% such that the frequency of recession in the simulation corresponds to the
% empirical frequency of recessions:
% Recession s_t < Rec_s_bar
rec_sim_ss = NaN(length(astsim_pf), 1);
for i = 1:length(astsim_pf)
    if astsim_pf(i) < Rec_s_bar
        rec_sim_ss(i) = 1;
    else
        rec_sim_ss(i) = 0;
    end
end