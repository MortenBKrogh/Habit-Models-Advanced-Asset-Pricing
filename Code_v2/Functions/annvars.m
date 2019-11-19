function [alndctsim astsim alnpctsim alnrtsim alnrfsim asdlnrtsim alnchpsim ...
    alnysim aelnrcbsim asdlnrcbsim atesterf aelnrtsim]=annvars(dc,lnpc,er,elnr,sdr,sdlnr,elnrcb,sdlnrcb,lny,lnrf1)
% Annualising and preparing data from the simulation "simvars.m". Returns
% various series of interest. Returns, PC and DC ratio, std, bond returns
% etc.

global tsc bondsel ann ncalc
% Simulating series
[stsim vtsim lndctsim lnpctsim lnrtsim lnrfsim ertsim elnrtsim sdrtsim...
    sdlnrtsim elnrcbsim sdlnrcbsim lnysim lnrcbsim testerf]=simvars(dc,lnpc,er,elnr,sdr,sdlnr,elnrcb,sdlnrcb,lny,lnrf1);

T = size(stsim,1);

%% Consumption
if ann == 1
    alndctsim=lndctsim;
else
    alnctsim = cumsum(lndctsim);    
    % Monthly logs
    alnctsim = log(chgfreq(exp(alnctsim),tsc,tsc,0));
    alndctsim = alnctsim(2:size(alnctsim,1))-alnctsim(1:(size(alnctsim,1)- 1));
    
end

%% s_t
if T > 1
    astsim = chgfreq(stsim(2:T),1,tsc,0);
    astsim = astsim(2:size(astsim,1)); 
end

%% P/C Ratio
if size(lnpctsim,1) > 1
    alnpctsim = chgfreq(lnpctsim(2:T),1,tsc,0)-log(tsc);
    alnpctsim = alnpctsim(2:size(alnpctsim,1));
end

%% Yearly Returns
if size(lnrtsim,1) > 1
    alnrtsim = chgfreq(lnrtsim,tsc,tsc,0);
    alnrtsim = alnrtsim(2:size(alnrtsim,1));
end

% Expected returns:
if size(elnrtsim,1) > 1
    aelnrtsim = chgfreq(elnrtsim,tsc,tsc,0);
    aelnrtsim = aelnrtsim(2:size(aelnrtsim,1));
end

% Risk free rate
if size(lnrfsim,1) > 1
    alnrfsim = chgfreq(lnrfsim(1:T-1),tsc,tsc,0);
    alnrfsim = alnrfsim(2:size(alnrfsim,1));
end
% Interpolation
if size(testerf,1) > 1
    atesterf = chgfreq(testerf(1:T-1),tsc,tsc,0);
    atesterf = atesterf(2:size(atesterf,1));
end
%% Conditional deviation of returns
if size(sdlnrtsim,1) > 1
    asdlnrtsim = chgfreq(sdlnrtsim,tsc,tsc,0);
    asdlnrtsim = asdlnrtsim(2:size(asdlnrtsim,1));
end
%% Price evolution
if size(lnpctsim,1) > 1
    lnchpsim = lnpctsim(2:T)-lnpctsim(1:T-1)+lndctsim;
    alnchpsim = chgfreq(lnchpsim,tsc,tsc,0);
    alnchpsim = alnchpsim(2:size(alnchpsim,1));
end
%% Yields
if size(lnysim,1) > 1
    for i=1:length(bondsel)+1
     alnysim(:,i) = chgfreq(lnysim(1:T-1,i),tsc,tsc,0);
    end
    alnysim = alnysim(2:size(alnysim,1),:); 
end
%% Bonds
% Mean returns
if size(elnrcbsim,1) > 1
    for i=1:length(bondsel)+1
       aelnrcbsim(:,i) = chgfreq(elnrcbsim(1:T-1,i),tsc,tsc,0);
    end
    
    aelnrcbsim = aelnrcbsim(2:size(aelnrcbsim,1),:);
end

% Deviations
if size(sdlnrcbsim,1) > 1
    for i=1:length(bondsel)
        asdlnrcbsim(:,i) = chgfreq(sdlnrcbsim(1:T-1,i+1),tsc,tsc,0);
    end
    
    asdlnrcbsim = asdlnrcbsim(2:size(asdlnrcbsim,1),:);
end
end