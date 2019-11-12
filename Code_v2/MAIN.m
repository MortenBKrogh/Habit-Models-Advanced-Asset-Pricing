clear all
format long
tic

addpath('Functions');
addpath('Data');
addpath('Workspaces');
addpath('Calibration');
%% 
% Calibration Choice
calib=2;    % 1 - Campbell & Cochrane (1999)
            % 2 - Krogh & Jensen (2019)
            
%%
global g sig delta phi gamma S_bar s_bar S_max s_max tsc sg B maxcb ncalc ...
    bondsel rhow seedval verd debug ann lnpca con sig_w lnpda PD_Claim

% Initialization

if calib == 1
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
    PD_Claim = 0; % 1 = PD_claim % 0 = PC_Claim
end

if calib == 2
    Pars = Krogh_Jensen_Calibration; % Change to Krogh_Jensen_Calibration also file name..
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
    PD_Claim = 0; % 1 = PD_claim % 0 = PC_Claim

end

PD_Claim_init = PD_Claim;
rho = (-1:.1:1);
S_bar=sig*sqrt(gamma/(1-phi-B/gamma));
s_bar = log(S_bar);
s_max = s_bar + (1-S_bar^2)/2;
S_max = exp(s_max);
delta=exp(gamma*g-.5*((1-phi)*gamma-B)-rf0); % Equation (12) in paper C&C- 1999.
szgrid=15;

ncalc = 100000; % Number if simulations
bondsel = [1 2 3 4 5 7 10 20]; % Maturity of bonds simulated
maxcb = max(bondsel);
seedval = 123;
chk = 1;

flag1 = 0; % Simulation flag
flag2 = 1; % 1 Simulation of yearly data, 0 of quarterly
con = 0;   % Interpolation

%% Grid def
sg = mkgrids(szgrid,0);
S=exp(sg);

%% PD- & PC-ratio

PD_Claim = 0; % 1 = PD_claim % 0 = PC_Claim
    
    figure;
    [lnpca ctrindx]=findlpc(sig,g,s_bar);
    PC_ratio=exp(lnpca);
    lnpca_pf=lnpca;
    plot(S,PC_ratio/tsc,'red'); % Annulized P/C-curve
    hold on;
    
PD_Claim = 1; % 1 = PD_claim % 0 = PC_Claim    

    [lnpda dtrindx]=findlpc(sig,g,s_bar);
    PD_ratio=exp(lnpda);
    plot(S,PD_ratio/tsc,'blue'); % Annulized P/C-curve
    legend('PC-Ratio', 'PD-Ratio')
    hold off;
    lnpda_pf=lnpda;
    
% reset PD_Claim to initial value we only changed it to make the plot
PD_Claim = PD_Claim_init;


%% Find expected returns and conditional deviations of consumption clain
verd=0;
% Fixed point method
    [er_pf elnr_pf sdr_pf sdlnr_pf lnrf_pf lnrf1_pf lny_pf elnrcb_pf sdlnrcb_pf slpmv_pf] = finders(sg);
    
%% Adjustments of inputs for simulation
if flag1 == 0 && flag2 == 0
    dc = 0;
elseif flag1 == 1 && flag2 == 0
    dc = exp(qlndc);
elseif flag1 == 0 && flag2 == 1
    dc = 0;
elseif flag1 == 1 && flag2 == 1
    dc = exp(alndc);
end

%% Simulation of time-series

randn('seed',seedval);
[alndctsim_pf astsim_pf alnpctsim_pf alnrtsim_pf alnrfsim_pf asdlnrtsim_pf ...
        alnchpsim_pf alnysim_pf aelnrcbsim_pf asdlnrcbsim_pf atesterfsim_pf] ...
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
    
    table = zeros(13,1);
    table(1,1) = Edc_pf;
    table(2,1) = Stdc_pf;
    table(3,1) = Erf_pf;
    table(4,1) = Stdrf_pf;
    table(5,1) = Shpr_pf;
    table(6,1) = ShpR_pf;
    table(7,1) = Eexrett_pf;
    table(8,1) = Stdexrett_pf;
    table(9,1) = Ep_d_pf;
    table(10,1) = Stdp_d_pf;
    table(11,1)= S_max;
    table(12,1)= S_bar;
    table(13,1)= delta^tsc;
%% SDF Simulation

randn('seed',seedval);
    [stsim, vtsim lndctsim lnpctsim lnrtsim lnrfsim ertsim elnrtsim sdrtsim...
        sdlnrtsim elnrcbsim sdlnrcbsim lnysim lnrcbsim testerfsim]=...
        simvars(dc,lnpca_pf,er_pf,elnr_pf,sdr_pf,sdlnr_pf,elnrcb_pf,sdlnrcb_pf,lny_pf ,lnrf1_pf);
    % Stochastic discount factor
    SDFus = delta*exp(-g*gamma)*exp(-gamma*vtsim).*exp(- gamma*(stsim(2:length(stsim))...
        -stsim(1:length(stsim)-1)));
    
    for k=1:length(rho)
        [stx vtx lndctx lnrfx]=simulacorr(rho(k));
        SDFx(:,k)= delta*exp(-g*gamma)*exp(-gamma*vtx).*exp(- gamma*(stx(2:length(stx))...
            -stx(1:length(stx)-1)));
        stsimx(:,k)=stx;
        vtsimx(:,k)=vtx;
        lndctsimx(:,k)=lndctx;
        lnrfsimx(:,k)=lnrfx;
    end
    %% 
    ts1 = struct();
    ts1.S_t           = astsim_pf;
    ts1.deltac       = alndctsim_pf;
    ts1.pcratio       = alnpctsim_pf;
    ts1.ExPostReturns = alnrtsim_pf;
    ts1.RiskFreeRate  = alnrfsim_pf;
    ts1.Prices        = alnchpsim_pf;
    ts1.stdReturns    = asdlnrtsim_pf;
    ts1 = struct2table(ts1);
    %writetable(ts1)
    %% Real FX-rates
    deltaq = log(SDFx) - kron(ones(1,size(SDFx,2)),log(SDFus));

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
rec_emp_percentage = sum(NBER_REC.data(idx_from:idx_to,1)) / length(NBER_REC.data(idx_from:idx_to,1))

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

    
%% Regression
   
    
    load gong
audioplayer(y,Fs);
play(ans)
toc

