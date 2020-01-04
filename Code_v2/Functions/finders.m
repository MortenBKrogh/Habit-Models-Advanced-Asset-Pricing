function [er elnr sdr sdlnr lnrf lnrf1 lny elnrcb sdlnrcb slpmv] = finders(sg)

%Procedure that calculates expected returns on consumer assets. elenos will 
% provide E (R), SD (R), lnrf, the mean boundary curvature of the variance 
% given by the variable (slpmv) and integrate the Sharpe ratio of the
% consumer asset vector.                                                   %
% ----------------------------------------------------------------------- %

global g gamma sig phi s maxcb s_bar delta tsc lnpcb matur PD_Claim

% Medium-Variance Boundary Slope                              %
%                                                                         %
% slpmv = (exp((gamma*sig)^2.*(1+lambda(sg)).^2)-1).^.5                   %
% ----------------------------------------------------------------------- %

slpmv = (exp((gamma*sig)^2*(1+lambda(sg)).^2)-1).^(0.5);

%% Term interest rate structure given by:                            %
%                                                                         %
% lnrf = -ln(delta) + gamma*g -                                           %
% gamma*(1-phi)*(s{t}-s_bar)-.5((gamma*sig)^2)*(1+lambda(s{t}))^2
% - B*(sg - s_bar) --> Vari?vel no estado
% ----------------------------------------------------------------------- %

lnrf = -log(delta) + gamma*g - gamma*(1-phi)*(sg-s_bar)...
    - 0.5*(gamma*sig*(1+lambda(sg))).^2;

%% Bonds
% Matrix of all bond prices. Your dimension will be N(sg) x (maxcb*tsc)
lnpcb = [];
lnpcb(:,1) = -lnrf;

lnp = zeros(size(sg,1),1);

for j = 2:(maxcb*tsc)
    for i = 1:length(sg)
        s = sg(i);
        lnp(i) = log(GaussLegendre(@intpcb,abs(sig)*(-8),abs(sig)*8,40));
    end
    lnpcb = cat(2,lnpcb,lnp);
end

% Yields
lny = - ...
    lnpcb./kron(ones(size(sg,1),1),linspace(1/tsc,(maxcb*tsc)/tsc,(maxcb*tsc)));

%% Expected Returns and Standard Deviations                                   %
% ----------------------------------------------------------------------- %

lnrf1 = zeros(size(sg,1),1);

er = zeros(size(sg,1),1);
elnr = zeros(size(sg,1),1);
sdr = zeros(size(sg,1),1);
sdlnr = zeros(size(sg,1),1);
elnrcb = zeros(size(sg,1),maxcb*tsc);       % zero-coupon bonds %
sdlnrcb = zeros(size(sg,1),maxcb*tsc);

for i=1:size(sg,1)
    s = sg(i);
    
    
    lnrf1(i)= - log(GaussLegendre(@intemrs,abs(sig)*(-8),abs(sig)*8,40));
    if PD_Claim == 0
    er(i)= GaussLegendre(@inter,abs(sig)*(-8),abs(sig)*8,40);
    sdr(i) = GaussLegendre(@inter2,abs(sig)*(-8),abs(sig)*8,40);
    else
    er(i) = GaussLegendre(@interd,abs(sig)*(-8),abs(sig)*8,40);
    sdr(i) = GaussLegendre(@inter2d,abs(sig)*(-8),abs(sig)*8,40);    
    end
    elnr(i)= GaussLegendre(@intelnr,abs(sig)*(-8),abs(sig)*8,40);
    sdr(i) = (sdr(i) - er(i).^2).^(.5);
    sdlnr(i) = GaussLegendre(@intelnr2,abs(sig)*(-8),abs(sig)*8,40);
    
    % Bonds
    matur = maxcb*tsc; elnrcb(i,1) = lnrf(i);
    
    while matur >= 2
        elnrcb(i,matur) = GaussLegendre(@intelnrcb,abs(sig)*(-8),abs(sig)*8, 40); 
        sdlnrcb(i,matur) = GaussLegendre(@intelnr2,abs(sig)*(-8),abs(sig)*8, 40); 
        sdlnrcb(i,matur) = (sdlnrcb(i,matur) - elnrcb(i,matur).^2).^(.5);
        matur = matur - 1;
    end
    
end
end


