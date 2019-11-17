function out = ercbin(v)
% 	Expected Returns of the Term Structure                                %
global s sg lnpcb matur
out = exp(interp(strans(s,v),sg,lnpcb(:,matur-1))) ./ ...
    exp(interp(s,sg,lnpcb(:,matur)));
end