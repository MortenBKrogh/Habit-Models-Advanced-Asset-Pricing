function out = intpcb(v)

% Function that provides the price of bonds for each maturity.            %

global s sg lnpcb sig

out = pdf('norm',v,0,sig).*mrsinsd(v).*...
    exp(interp(strans(s,v),sg,lnpcb(:,size(lnpcb,2))))';

end