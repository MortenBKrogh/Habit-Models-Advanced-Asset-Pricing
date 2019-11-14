function [stsim vtsim lndctsim lnrfsim]=simulacorr(rho)

global ncalc gamma sig g phi delta B s_bar seedval

%% Simulation of shocks from the correlationcoefficient rho

T=ncalc;
randn('seed',seedval);
x = sig*randn(T,1);
y = sig*randn(T,1);

vtsim = rho*x + sqrt(1-rho^2)*y; 
lndctsim = g + vtsim;

%% Simulation log state (Surplus consumption ratio)

stsim = zeros(T+1,1);

stsim(1) = s_bar;        

for i=2:T+1
    
    stsim(i) = strans(stsim(i-1),vtsim(i-1));

end

%% Time variying log RF-rate

lnrfsim = -log(delta)+gamma*g-(gamma*(1-phi)-B)/2-B .*(stsim-s_bar); 

end