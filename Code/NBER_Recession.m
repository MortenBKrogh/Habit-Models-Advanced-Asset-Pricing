clear; clc;
% Load NBER Recession data from 1854-12-01 to 2019-10-01
% The USREC.csv is monthly observed.
% For updated data see
NBER_REC = importdata('USREC.csv');

% Define period yyyy-mm-dd
from = '1950-01-01';
to   = '2018-12-01';

% find indexes
idx_from = find(NBER_REC.textdata(:,1)==string(from)) - 1;
idx_to   = find(NBER_REC.textdata(:,1)==string(to)) - 1;

% Calculate percentage of the time the economy is in recession
rec_emp_percentage = sum(NBER_REC.data(idx_from:idx_to,1)) / length(NBER_REC.data(idx_from:idx_to,1))

%% Now we look at how to define the recession dummey to match the empirical findings
% Load the saved workspace
load('Reduced_workspace.mat');
clearvars -except astsim_pf s_bar from to idx_from idx_to NBER_REC rec_emp_percentage 

% define recession dummey as when s_t < s_bar
rec_sim_ss = NaN(length(astsim_pf), 1);

for i = 1:length(astsim_pf)
   
    if astsim_pf(i) < s_bar
        rec_sim_ss(i) = 1;
    else 
        rec_sim_ss(i) = 0;
        
    end
    
end
rec_sim_ss_percentage = sum(rec_sim_ss(:)==1) / length(rec_sim_ss)

% We see a big difference in out simulated economy compared to what
% empirical has been the case, thus we want to find a threshold which is
% lower than steady state value of surplus consumption growth. 

rec_sim_2 = NaN(length(astsim_pf), 1);

%%
vec = s_emp_recession(s_bar, astsim_pf)

p_rec = percentage(vec)

% finding value of s_bar that makes the amount of recessions mathc the
% empirical amount of recessions. 

res = @(s_bar, s, emp_rec) percentage(s_emp_recession(s_bar, s)) - emp_rec
%%
fplot(@(s_bar) percentage(s_emp_recession(s_bar, astsim_pf)))
%%
fplot(@(s_bar) res(s_bar, astsim_pf, rec_emp_percentage), [-15, 10])