function [stsim vtsim lndctsim lnpctsim lnrtsim lnrfsim ertsim elnrtsim sdrtsim...
    sdlnrtsim elnrcbsim sdlnrcbsim lnysim lnrcbsim testerf]=simvars(dc,lnpca,er,elnr,sdr,sdlnr,elnrcb,sdlnrcb,lny,lnrf1)

%
% This routine simulates the most important time-series of this model     %
% from a chosen calibration                                               %
% Simulating:                                                             %
% - s=log(S);                                                             %
% - P/C;                                                                  %
% - R{t+1};                                                               %
% - E{t}[R{t+1}];                                                         %
% - SD{t}[R{t+1}];                                                        %
% - Rf{t+1};                                                              %
% - corr(Rf{t+1},Cons_t)                                                  %
% - Bonds                                                                 %
%                                                                         %

global ncalc gamma sig sig_w g phi delta s_max s_bar sg maxcb tsc bondsel...
    PD_Claim rhow g

%% initialization
if dc == 0 
    T=ncalc;
    vtsim = sig*randn(T,1);
    wtsim = rhow * sig_w / sig * vtsim + sig_w * (1 - rhow ^2) ^ 0.5 * randn(T,1);
    if PD_Claim == 0
        lndctsim = g + vtsim;
    else
        lndctsim = g + wtsim;
    end
else
    if min(dc) <= 0
        
        disp ('simvars: You entered the consumption growth log.');
        disp ('You need to enter consumption growth data');
        disp ('gross, ie neither log nor net growth.');
        
    end
    T = length(dc);
    lndctsim = log(dc); 
end

%% Simulation of log(S_t)

stsim = zeros(T+1,1);

stsim(1) = s_bar;           % Starting the process at SS

if PD_Claim == 0
    for i=2:T+1
        if strans(stsim(i-1),vtsim(i-1)) <= s_max
            stsim(i) = strans(stsim(i-1),vtsim(i-1));
        else
            stsim(i)=(1-phi)*s_bar+phi*stsim(i-1);
        end
    end
else
    for i=2:T+1
        if strans(stsim(i-1),wtsim(i-1)) <= s_max
            stsim(i) = strans(stsim(i-1),vtsim(i-1));
        else
            stsim(i)=(1-phi)*s_bar+phi*stsim(i-1);
        end
    end
end
%% PC ratio                                           % 
lnpctsim = interp(stsim,sg,lnpca)';

%% ex-post Returns                                                        %
%                                                                         % 
%                        R = (C'/C){(1+(P/C)')/(P/C)}                     % 
% ----------------------------------------------------------------------- %
if PD_Claim == 0
lnrtsim = log(1+exp(lnpctsim(2:T+1))) - lnpctsim(1:T) + lndctsim;
else
lnrtsim =log( 1+exp(lnpctsim(2:T+1)) ) - lnpctsim(1:T) + lndctsim;
end
    %%
%% potential time varying RF-rate

lnrfsim = -log(delta) + gamma*g - gamma*(1-phi)*(stsim-s_bar)... 
    - 0.5*(gamma*sig*(1+lambda(stsim))).^2;

%% Expected Returns and Expected deviations 

testerf = interp(stsim,sg,lnrf1)';
ertsim = interp(stsim,sg,er)';
elnrtsim = interp(stsim,sg,elnr)'; 
sdrtsim = interp(stsim,sg,sdr)'; 
sdlnrtsim = interp(stsim,sg,sdlnr)';

%% Treasury Bill return 90 days

elnrcbsim = interp(stsim,sg,elnrcb(:,1))'; % Expected Returns on Tbill
sdlnrcbsim = interp(stsim,sg,sdlnrcb(:,1))'; 
lnysim = interp(stsim,sg,lny(:,1))';       % Bond yields
lny2sim = zeros(size(lnysim,1),1);

for i = 2:(maxcb*tsc)
    if find(i == bondsel*tsc)
        
        lnysim = cat(2,lnysim,interp(stsim,sg,lny(:,i))');
        lny2sim = cat(2,lny2sim,interp(stsim,sg,lny(:,i-1))'); 
        elnrcbsim = cat(2,elnrcbsim,interp(stsim,sg,elnrcb(:,i))'); 
        sdlnrcbsim = cat(2,sdlnrcbsim,interp(stsim,sg,sdlnrcb(:,i))');
    
    end
end

% Returns on bonds with maturities = [1 2 4 8 12 16 20] 

lnrcbsim = cat(1,0,lnysim(1:T-1,1));
for i = 2:length(bondsel)
    lnrcbsim = cat(2,lnrcbsim,cat(1,0,(-lny2sim(2:T,i)*((bondsel(i)- 1/tsc))+...
        lnysim(1:T-1,i)*bondsel(i))/tsc));
end
end


