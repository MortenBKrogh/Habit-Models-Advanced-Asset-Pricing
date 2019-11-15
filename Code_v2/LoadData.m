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
load('Workspaces/Calibration','s_bar');
load('Workspaces/CC_PC_Claim_workspace','rec_emp_percentage')
%% Matching the empirical density
%Rec_s_bar = fzero(@(x) (integral(@q_s,-Inf,x) - rec_emp_percentage), s_bar-0.9);
Rec_s_bar = -2.75;
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