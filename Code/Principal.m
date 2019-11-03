clear all
format long
tic
%% 
% Choice of method
method = 0; % Method: 0 Fixed point
            % Method: 1 Wachter 2005
            % Method: 2 comparative
% Calibration Choice
calib=5;    % 0 - b>0
            % 1 - b<0
            % 2 - Campbell Cochrane (1999)
            % 3 - Paper Wachter (2005)
            % 4 - Working Paper Verdelhan (2008)
            % 5 - Krogh Jensen (2019)
            
%%
global g sig delta phi gamma S_bar s_bar S_max s_max tsc sg B maxcb ncalc ...
    bondsel rho seedval verd debug ann lnpca con
% Initialization
if calib == 0
    tsc = 4; % Interval 4 = quarter
    g=0.0228/tsc;
    sig=0.009/sqrt(tsc);
    rf0=0.0127/tsc;
    B=0.009;
    gamma=2.5;
    phi=0.902^(1/tsc);
    verd=0;
    ann=0;
end

if calib == 1
    tsc = 4;
    g=0.0219/tsc;
    sig=0.0202/sqrt(tsc);
    rf0=0.0098/tsc;
    phi=0.931^(1/tsc);
    gamma=2;
    B=-0.01;
    verd=1;
    ann=0;
end

if calib == 2
    tsc = 12;
    g=0.0189/tsc;
    sig=0.015/sqrt(tsc);
    rf0=0.0094/tsc;
    phi=0.87^(1/tsc);
    gamma=2;
    B=0;
    verd=0;
    ann=0;
end

if calib == 3
    tsc = 4;
    g=0.022/tsc;
    sig=0.0086/sqrt(tsc);
    rf0=0.0147/tsc;
    phi=(0.89)^(1/tsc);
    gamma=2;
    B=0.011;
    verd=0;
    ann=0;
end

if calib == 4
    tsc = 4;
    g=0.02115/tsc;
    sig=0.0102/sqrt(tsc);
    rf0=0.0136/tsc;
    phi=0.97^(1/tsc);
    gamma=2.4;
    B=-0.01;
    verd=1;
    ann=1;
end

if calib == 5
    Pars = Model_Calibration;
    tsc = 12;
    g = Pars.g/tsc;
    sig = Pars.sigma/sqrt(tsc);
    rf0 = Pars.rf/tsc;
    phi=Pars.Phi^(1/tsc);
    gamma=2;
    B=0;
    verd=0;
    ann=0;
end

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

if method == 0 || method == 2
    [lnpca ctrindx]=findlpc(sig,g,s_bar);
    PC_ratio=exp(lnpca);
    plot(S,PC_ratio/tsc,'k'); % Annulized P/C-curve
end
lnpca_pf=lnpca;

%% Find Value of P/C, serial method
if method == 1 || method == 2
    [W_PC_ratio]=WfindFn(sig,sg); 
    plot(S,W_PC_ratio/tsc,'k');
    lnpca_s=log(W_PC_ratio);
    lnpca=lnpca_s;
end
% Comparative method
if method == 2
    plot(S,W_PC_ratio/tsc,'r',S,PC_ratio/tsc,'g'); legend('Series method','Fixed-point method',2); xlabel('Consumption surplus ratio (S{t})'); ylabel('Price-consumption ratio (P{t}/C{t})'); comp=max(abs((W_PC_ratio-PC_ratio)./PC_ratio));
end

%% Find expected returns and conditional deviations of consumption clain
verd=0;
% Fixed point method
if method == 0 || method == 2
    [er_pf elnr_pf sdr_pf sdlnr_pf lnrf_pf lnrf1_pf lny_pf elnrcb_pf sdlnrcb_pf slpmv_pf] = finders(sg);
end
% Serial method
if method == 1 || method == 2
    [er_s elnr_s sdr_s sdlnr_s lnrf_s lnrf1_s lny_s elnrcb_s sdlnrcb_s slpmv_s] = finders(sg);
end
%%
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
if method == 0 || method == 2
    randn('seed',seedval);
    [alndctsim_pf astsim_pf alnpctsim_pf alnrtsim_pf alnrfsim_pf asdlnrtsim_pf ...
        alnchpsim_pf alnysim_pf aelnrcbsim_pf asdlnrcbsim_pf atesterfsim_pf]=annvars(dc,lnpca_pf,er_pf,elnr_pf,sdr_pf,sdlnr_pf,elnrcb_pf,sdlnrcb_pf,lny_pf,lnrf1_pf);
end
if method == 1 || method == 2
    randn('seed',seedval);
    [alndctsim_s astsim_s alnpctsim_s alnrtsim_s alnrfsim_s sdlnrtsim_s alnchpsim_s ...
        alnysim_s aelnrcbsim_s asdlnrcbsim_s atesterfsim_s]=annvars(dc,lnpca_s,er_s,elnr_s,sdr_s,...
        sdlnr_s,elnrcb_s,sdlnrcb_s,lny_s,lnrf1_s);
end

%% Statistics of interest
if method == 0 || method == 2
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
    Ep_d_pf = mean(alnpctsim_pf); % Log price consumption
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
end
if method == 1 || method == 2
    if ann == 1
        Edc_s = tsc*mean(alndctsim_s);
        Stdc_s = sqrt(tsc)*std(alndctsim_s);
    else
        Edc_s = mean(alndctsim_s);
        Stdc_s = std(alndctsim_s);
    end
    Erf_s = mean(alnrfsim_s);
    Stdrf_s = std(alnrfsim_s);
    
    Erfinterp_s = mean(atesterfsim_s);
    Stdrfinterp_s = std(atesterfsim_s);
    
    exrett_s = alnrtsim_s - alnrfsim_s; % Excess returns
    exrettinterp_s = alnrtsim_s - atesterfsim_s;
    Shpr_s = mean(exrett_s)/std(exrett_s); % Sharpe Ratio
    ShpR_s = mean(exp(alnrtsim_s)-exp(alnrfsim_s))/std(exp(alnrtsim_s)- exp(alnrfsim_s));
    Shprinterp_s = mean(exrettinterp_s)/std(exrettinterp_s);
    ShpRinterp_s = mean(exp(alnrtsim_s)-exp(atesterfsim_s))/std(exp(alnrtsim_s)- exp(atesterfsim_s));
    
    Eexrett_s = mean(exrett_s); % Mean excess returns
    Stdexrett_s = std(exrett_s); 
    Eexrettinterp_s = mean(exrettinterp_s);
    Stdexrettinterp_s = std(exrettinterp_s);
    
    Ep_d_s = mean(alnpctsim_s); % Mean log consumption
    Stdp_d_s = std(alnpctsim_s); 
    
    table = zeros(13,1);
    table(1,1) = Edc_s;
    table(2,1) = Stdc_s;
    table(3,1) = Erf_s;
    table(4,1) = Stdrf_s;
    table(5,1) = Shpr_s;
    table(6,1) = ShpR_s;
    table(7,1) = Eexrett_s;
    table(8,1) = Stdexrett_s;
    table(9,1) = Ep_d_s;
    table(10,1) = Stdp_d_s;
    table(11,1)= S_max;
    table(12,1)= S_bar;
    table(13,1)= delta^tsc;
end

if method == 2
    table(1,1) = Edc_pf; table(1,2)=Edc_s;
    table(2,1) = Stdc_pf; table(2,2) = Stdc_s;
    table(3,1) = Erf_pf; table(3,2) = Erf_s;
    table(4,1) = Stdrf_pf; table(4,2) = Stdrf_s;
    table(5,1) = Shpr_pf; table(5,2) = Shpr_s;
    table(6,1) = ShpR_pf; table(6,2) = ShpR_s;
    table(7,1) = Eexrett_pf; table(7,2) = Eexrett_s;
    table(8,1) = Stdexrett_pf; table(8,2) = Stdexrett_s;
    table(9,1) = Ep_d_pf; table(9,2) = Ep_d_s;
    table(10,1) = Stdp_d_pf; table(10,2) = Stdp_d_s;
    table(11,1)= S_max; table(11,2)= S_max;
    table(12,1)= S_bar; table(12,2)= S_bar;
    table(13,1)= delta^tsc; table(13,2)= delta^tsc;
end

%% SDF Simulation

if method == 0 || method == 2
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
    writetable(ts1)
    %% Real FX-rates
    deltaq = log(SDFx) - kron(ones(1,size(SDFx,2)),log(SDFus));
    %% Regression
    if con == 1
        for k=1:length(rho)
            diff(:,k)=lnrfsim-lnrfsimx(:,k);
            regressor=cat(2,ones(length(diff(:,k)),1),diff(:,k));
            [betas interv]=regress(deltaq(:,k),regressor(1:length(regressor)-1,:));
            Betas_pf(k,1)=rho(k);
            Betas_pf(k,2)=betas(1,1);
            Betas_pf(k,3)=betas(2,1);
            if k==1
                Interv_pf=interv;
            else
                Interv_pf=cat(1,Interv_pf,interv);
            end
        end
    elseif con == 0
        for k=1:length(rho)
            diff(:,k)=testerfsim-lnrfsimx(:,k);
            regressor=cat(2,ones(length(diff(:,k)),1),diff(:,k));
            [betas interv]=regress(deltaq(:,k),regressor(1:length(regressor)-1,:));
            Betas_pf(k,1)=rho(k);
            Betas_pf(k,2)=betas(1,1);
            Betas_pf(k,3)=betas(2,1);
            if k==1
                Interv_pf=interv;
            else
                Interv_pf=cat(1,Interv_pf,interv);
            end
        end
    end
end
%%
if method == 1 || method == 2
    randn('seed',seedval); 
    [stsim vtsim lndctsim lnpctsim lnrtsim lnrfsim ertsim elnrtsim sdrtsim...
        sdlnrtsim elnrcbsim sdlnrcbsim lnysim lnrcbsim testerfsim]=...
        simvars(dc,lnpca_s,er_s,elnr_s,sdr_s,sdlnr_s,elnrcb_s,sdlnrcb_s,lny_s,lnrf1_s );
    
    SDFus = delta*exp(-g*gamma)*exp(-gamma*vtsim).*exp(- gamma*(stsim(2:length(stsim))...
        - stsim(1:length(stsim)-1)));
    
    for k=1:length(rho)
        [stx vtx lndctx lnrfx]=simulacorr(rho(k),lnrf1_s);
        SDFx(:,k)= delta*exp(-g*gamma)*exp(-gamma*vtx).*exp(- gamma*(stx(2:length(stx))...
            -stx(1:length(stx)-1)));
        stsimx(:,k)=stx;
        vtsimx(:,k)=vtx;
        lndctsimx(:,k)=lndctx;
        lnrfsimx(:,k)=lnrfx;
    end
    deltaq = log(SDFx) - kron(ones(1,size(SDFx,2)),log(SDFus));
    if con == 1
        for k=1:length(rho)
            diff(:,k)=lnrfsim-lnrfsimx(:,k);
            regressor=cat(2,ones(length(diff(:,k)),1),diff(:,k));
            [betas interv]=regress(deltaq(:,k),regressor(1:length(regressor)-1,:));
            Betas_s(k,1)=rho(k);
            Betas_s(k,2)=betas(1,1);
            Betas_s(k,3)=betas(2,1);
            if k==1
                Interv_s=interv;
            else
                Interv_s=cat(1,Interv_s,interv);
            end
        end
    elseif con == 0
        for k=1:length(rho)
            diff(:,k)=testerfsim-lnrfsimx(:,k);
            regressor=cat(2,ones(length(diff(:,k)),1),diff(:,k));
            [betas interv]=regress(deltaq(:,k),regressor(1:length(regressor)-1,:));
            Betas_s(k,1)=rho(k);
            Betas_s(k,2)=betas(1,1);
            Betas_s(k,3)=betas(2,1);
            if k==1
                Interv_s=interv;
            else
                Interv_s=cat(1,Interv_s,interv);
            end
        end
    end
end
if method == 2
    Betas = cat(2,Betas_pf,NaN*ones(size(Betas_pf,1),1));
    Betas = cat(2,Betas,Betas_s);
    Intervalo = cat(2,Interv_pf,NaN*ones(size(Interv_pf,1),1));
    Intervalo = cat(2,Intervalo,Interv_s);
end
%{
bn=zeros(length(bondsel)-1,2);
bnconf=zeros(length(bondsel)-1,4);
for i=2:length(bondsel)
    chgyield=lnysim(2:size(lnysim,1),i-1)-lnysim(1:size(lnysim,1)-1,i);
    spreadyield=(lnysim(1:size(lnysim,1)-1,i)-1 - nysim(1:size(lnysim,1)- 1,1))./(bondsel(i)-1);
    
    [aux1 aux2]=regress(chgyield,[ones(length(lnysim)-1,1) spreadyield]);
    bn(i-1,:)=aux1';
    bnconf(i-1,:)=reshape(aux2,1,4);
end
%}
load gong
audioplayer(y,Fs);
play(ans)
toc