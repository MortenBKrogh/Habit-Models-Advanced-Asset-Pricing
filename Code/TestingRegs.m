clear
clc
%%
Dat = readtable('ts1.txt');

s_t = Dat.S_t; % Is logged
pc_t = Dat.PCratio; % Is logged
ret_t = Dat.ExPostReturns; % Is logged i think
rfrate = Dat.RiskFreeRate; %
Net_ret_t = ret_t - rfrate;
%
ret_t = Net_ret_t;
%% training:
S_t_train = s_t(1:5000,1);
pc_t_train = pc_t(1:5000, 1);
ret_t_train = ret_t(2:5001, 1);
%%
pred = [ones(size(pc_t_train, 1), 1), pc_t_train, S_t_train];
[beta, std,~,~,stats] = regress(ret_t_train, pred);
Variables = {'Constant' 'pc_t-1' 's_t-1'};
Vect = [Variables' num2cell(beta) num2cell(std(:,1)) num2cell(std(:,2))];
% Format output
fprintf('    -------------------------------------| Model |-------------------------------------\n')
fprintf('    ###################################################################################\n')
fprintf('    ##  r_t = constant + beta_1 * log(P/C_{t-1}) + beta_2 * log(S_{t-1}) + epsilon_t ##\n')
fprintf('    ###################################################################################\n')
fprintf('    \n')
fprintf('    ----------------| Model Effects |----------------\n')
fprintf('    Variable      Coef.        Std_L.       Std_U  \n')
disp(Vect)
fprintf('    \n')
fprintf('    --------| Model Diagnostics |-------\n')
fprintf('    R2      F-stat      P-val     ResVar \n')
disp(stats)
%%
fprintf('  \n')
fprintf('  \n')
fprintf('  \n')
pred = [ones(size(pc_t_train, 1), 1), pc_t_train];
[beta, std,~,~,stats] = regress(ret_t_train, pred);
Variables = {'Constant' 'pc_t-1'};
Vect = [Variables' num2cell(beta) num2cell(std(:,1)) num2cell(std(:,2))];
% Format output
fprintf('    --------------------------| Model |------------------------\n')
fprintf('    ###########################################################\n')
fprintf('    ##  r_t = constant + beta_1 * log(deltaC_{t-1}) + epsilon_t ##\n')
fprintf('    ###########################################################\n')
fprintf('    \n')
fprintf('    ----------------| Model Effects |----------------\n')
fprintf('    Variable      Coef.        Std_L.       Std_U  \n')
disp(Vect)
fprintf('    \n')
fprintf('    --------| Model Diagnostics |-------\n')
fprintf('    R2      F-stat      P-val     ResVar \n')
disp(stats)