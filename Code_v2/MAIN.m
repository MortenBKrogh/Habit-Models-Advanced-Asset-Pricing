clear all
clc
format long
tic
addpath('Functions');
addpath('Data');
addpath('Workspaces');
addpath('Calibration');
addpath('Figures');
addpath('Tables');
%% Defining Globals
global g sig delta phi gamma S_bar s_bar S_max s_max tsc sg B maxcb ncalc ...
    bondsel rhow seedval verd debug ann lnpca con sig_w lnpda PD_Claim Regressions

%% Choices for solution methods
% Calibration Choice
calib=1;           % 0 - Campbell & Cochrane (1999)
                   % 1 - Krogh & Jensen (2019)

% Solution method:
PD_Claim = 1;      % 0 = Consumption Claim
                   % 1 = Dividend Claim
% Plots
Plots = 0;         % 0 = off
                   % 1 = on
% Update tables
Tables = 0;        % 0 = off
                   % 1 = on
% Regressions
Regressions = 0;   % 0 = off
                   % 1 = on
%% Initialization
if calib == 0
    tsc = 12;
    g=0.0189/tsc;
    sig=0.015/sqrt(tsc);
    rf0=0.0094/tsc;
    phi=0.87^(1/tsc);
    gamma=2;
    B=0;
    verd=0;
    ann=0;
    sig_w = 0.112/sqrt(tsc);
    rhow = 0.2;
end

if calib == 1
    Pars = Model_Calibration;
    tsc = 12;
    g = Pars.g/tsc;
    sig = Pars.sigma/sqrt(tsc);
    rf0 = Pars.rf/tsc;
    phi=Pars.Phi^(1/tsc);
    gamma=2;
    rhow = 0.2;
    B=0;
    verd=0;
    ann=0;
    sig_w = Pars.sigma_w/sqrt(tsc);
end

PD_Claim_init = PD_Claim;
rho = (-1:.1:1);

S_bar=sig*sqrt(gamma/(1-phi-B/gamma));
s_bar = log(S_bar);
s_max = s_bar + (1-S_bar^2)/2;
S_max = exp(s_max);
delta=exp(gamma*g-.5*((1-phi)*gamma-B)-rf0); % Equation (12) in paper C&C- 1999.

szgrid=15;

ncalc = 100000;                % Number of simulations
bondsel = [1 2 3 4 5 7 10 20]; % Maturity of bonds simulated
maxcb = max(bondsel);
seedval = 123;

chk = 1;
flag2 = 1; % 1 Simulation of yearly data, 0 of quarterly
con = 0;   % Interpolation

%% Grid def
sg = mkgrids(szgrid);
S=exp(sg);

%% PD- & PC-ratio
if Plots
 PD_Claim = 0;
 [output_lnpca ctrindx]=findlpc(sig,g,s_bar);
 PC_ratio=exp(output_lnpca);
 lnpca_pf=output_lnpca;

 PD_Claim = 1;
 [output_lnpda dtrindx]=findlpc(sig,g,s_bar);
 PD_ratio=exp(output_lnpda);
 lnpda_pf=output_lnpda;
end
clear lnpca lnpc
global lnpca lnpc
% reset PD_Claim to initial value we only changed it to make the plot
 PD_Claim = PD_Claim_init;

if PD_Claim == 0
    [output_lnpca ctrindx]=findlpc(sig,g,s_bar);
    lnpca = output_lnpca;
    lnpca_pf = output_lnpca;
else 
    [output_lnpda dtrindx]=findlpc(sig,g,s_bar);
    lnpca = output_lnpda;
    lnpca_pf = output_lnpda;
end
%% Find expected returns and conditional deviations of consumption claim
verd=0;
% Fixed point method
[er_pf elnr_pf sdr_pf sdlnr_pf lnrf_pf lnrf1_pf lny_pf elnrcb_pf sdlnrcb_pf slpmv_pf] = finders(sg);

%% Adjustments of inputs for simulation
dc = 0;
%% Simulation of time-series
[alndctsim_pf astsim_pf alnpctsim_pf alnrtsim_pf alnrfsim_pf asdlnrtsim_pf ...
    alnchpsim_pf alnysim_pf aelnrcbsim_pf asdlnrcbsim_pf atesterfsim_pf aelnrtsim] ...
    =annvars(dc,lnpca_pf,er_pf,elnr_pf,sdr_pf,sdlnr_pf,elnrcb_pf,sdlnrcb_pf,lny_pf,lnrf1_pf);
%% Statistics of interest


if ann == 1
    Edc_pf = tsc*mean(alndctsim_pf);
    Stdc_pf = sqrt(tsc)*std(alndctsim_pf);
else
    Edc_pf = mean(alndctsim_pf);
    Stdc_pf = std(alndctsim_pf);
end

Erf_pf = mean(alnrfsim_pf); % mean log riskfree rate
Stdrf_pf = std(alnrfsim_pf); % sd log RF-rate
Erfinterp_pf = mean(atesterfsim_pf);
Stdrfinterp_pf = std(atesterfsim_pf);

exrett_pf = alnrtsim_pf - alnrfsim_pf; % Excess returns
exrettinterp_pf = alnrtsim_pf - atesterfsim_pf;

Shpr_pf = mean(exrett_pf)/std(exrett_pf); % Sharpe ratio of log returns
ShpR_pf = mean(exp(alnrtsim_pf)-exp(alnrfsim_pf))/std(exp(alnrtsim_pf)- exp(alnrfsim_pf));
Shprinterp_pf = mean(exrettinterp_pf)/std(exrettinterp_pf);
ShpRinterp_pf = mean(exp(alnrtsim_pf)- exp(atesterfsim_pf))/std(exp(alnrtsim_pf)-exp(atesterfsim_pf));
Eexrett_pf = mean(exrett_pf); % Mean excess log returns
Stdexrett_pf = std(exrett_pf); % SD excess log returns)
Eexrettinterp_pf = mean(exrettinterp_pf);
Stdexrettinterp_pf = std(exrettinterp_pf);
Ep_d_pf = mean(alnpctsim_pf); % Log price/consumption
Stdp_d_pf = std(alnpctsim_pf);

if PD_Claim == 0
    PC_Claim_Sim_mom = struct();
    PC_Claim_Sim_mom.MeanConsGrowth       = Edc_pf;
    PC_Claim_Sim_mom.StdConsGrowth        = Stdc_pf;
    PC_Claim_Sim_mom.MeanRiskFreeRate     = Erfinterp_pf;
    PC_Claim_Sim_mom.StdRiskFreeRate      = Stdrfinterp_pf;
    PC_Claim_Sim_mom.logSharperatio       = Shprinterp_pf;
    PC_Claim_Sim_mom.Sharperatio          = ShpRinterp_pf;
    PC_Claim_Sim_mom.MeanExcessReturns    = Eexrettinterp_pf;
    PC_Claim_Sim_mom.StdExcessReturns     = Stdexrettinterp_pf;
    PC_Claim_Sim_mom.MeanPriceDividend    = Ep_d_pf;
    PC_Claim_Sim_mom.StdPriceDividend     = Stdp_d_pf;
    PC_Claim_Sim_mom.S_max                = S_max;
    PC_Claim_Sim_mom.S_bar                = S_bar;
    PC_Claim_Sim_mom.delta                = delta^tsc;
    PC_Claim_Sim_mom = struct2table(PC_Claim_Sim_mom);
    writetable(PC_Claim_Sim_mom)
elseif PD_Claim == 1
    PD_Claim_Sim_mom = struct();
    PD_Claim_Sim_mom.MeanDivGrowth        = Edc_pf;
    PD_Claim_Sim_mom.StdDivGrowth         = Stdc_pf;
    PD_Claim_Sim_mom.MeanRiskFreeRate     = Erfinterp_pf;
    PD_Claim_Sim_mom.StdRiskFreeRate      = Stdrfinterp_pf;
    PD_Claim_Sim_mom.logSharperatio       = Shprinterp_pf;
    PD_Claim_Sim_mom.Sharperatio          = ShpRinterp_pf;
    PD_Claim_Sim_mom.MeanExcessReturns    = Eexrettinterp_pf;
    PD_Claim_Sim_mom.StdExcessReturns     = Stdexrettinterp_pf;
    PD_Claim_Sim_mom.MeanPriceDividend    = Ep_d_pf;
    PD_Claim_Sim_mom.StdPriceDividend     = Stdp_d_pf;
    PD_Claim_Sim_mom.S_max                = S_max;
    PD_Claim_Sim_mom.S_bar                = S_bar;
    PD_Claim_Sim_mom.delta                = delta^tsc;
    PD_Claim_Sim_mom = struct2table(PD_Claim_Sim_mom);
    writetable(PD_Claim_Sim_mom)
end
%% SDF Simulation
rng(24,'twister')
[stsim, vtsim lndctsim lnpctsim lnrtsim lnrfsim ertsim elnrtsim sdrtsim...
    sdlnrtsim elnrcbsim sdlnrcbsim lnysim lnrcbsim testerfsim]=...
    simvars(dc,lnpca_pf,er_pf,elnr_pf,sdr_pf,sdlnr_pf,elnrcb_pf,sdlnrcb_pf,lny_pf ,lnrf1_pf);
% Stochastic discount factor
SDFus = delta*exp(-g*gamma)*exp(-gamma*vtsim).*exp(- gamma*(stsim(2:length(stsim))...
    -stsim(1:length(stsim)-1)));
%% Data Table
if PD_Claim == 0
    PC_Claim_Sim_dat = struct();
    PC_Claim_Sim_dat.S_t           = astsim_pf;
    PC_Claim_Sim_dat.deltac        = alndctsim_pf;
    PC_Claim_Sim_dat.pcratio       = alnpctsim_pf;
    PC_Claim_Sim_dat.ExPostReturns = alnrtsim_pf;
    PC_Claim_Sim_dat.RiskFreeRate  = alnrfsim_pf;
    PC_Claim_Sim_dat.Prices        = alnchpsim_pf;
    PC_Claim_Sim_dat.stdReturns    = asdlnrtsim_pf;
    PC_Claim_Sim_dat = struct2table(PC_Claim_Sim_dat);
    writetable(PC_Claim_Sim_dat)
elseif PD_Claim == 1
    PD_Claim_Sim_dat = struct();
    PD_Claim_Sim_dat.S_t           = astsim_pf;
    PD_Claim_Sim_dat.deltac        = alndctsim_pf;
    PD_Claim_Sim_dat.pcratio       = alnpctsim_pf;
    PD_Claim_Sim_dat.ExPostReturns = alnrtsim_pf;
    PD_Claim_Sim_dat.RiskFreeRate  = alnrfsim_pf;
    PD_Claim_Sim_dat.Prices        = alnchpsim_pf;
    PD_Claim_Sim_dat.stdReturns    = asdlnrtsim_pf;
    PD_Claim_Sim_dat = struct2table(PD_Claim_Sim_dat);
    writetable(PD_Claim_Sim_dat)
end
%% Construction of Indicator of recession
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
%% Regression
if Regressions == 1
returns = alnrtsim_pf;     % Returns
PD_regress = alnpctsim_pf; % PD / PC

h   = 0;                        % Forecast Horizon 0 = in-sample regression

y   = returns(1+h:end,1);
x   = [ones(length(returns(1:end-h,:)), 1),  ...         % vector of Ones
    rec_sim_ss(1:end-h,:) .* PD_regress(1:end-h,1), ...  %    I_rec_t *PD_t
    (1-rec_sim_ss(1:end-h,:)) .* PD_regress(1:end-h,1)]; % (1-I_rec_t)*PD_t
reg = nwest(y,x,0);
end
%% Finish
if Plots == 1
    Figures_CC1998;
end

if Tables == 1
    Table_Generator;
end

if calib == 1
    if PD_Claim == 0
        save('Workspaces/PC_Claim_workspace');
    else
        save('Workspaces/PD_Claim_workspace');
    end
else
    if PD_Claim == 0
        save('Workspaces/CC_PC_Claim_workspace');
    else
        save('Workspaces/CC_PD_Claim_workspace');
    end
end

%%
%load gong
%audioplayer(y,Fs);
%play(ans)
mean(aelnrtsim)-Erf_pf
toc
