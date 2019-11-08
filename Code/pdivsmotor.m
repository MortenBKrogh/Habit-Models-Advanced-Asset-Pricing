function [inside] = pdivsmotor(w)
% Procedure to be used for numeric integration when calculating
% the fixed point. It has as argument only w ~ N (0, sig_w). Returns VALUE of
% P / D for each iteration over the current value of s {t} in each iteration.
% ----------------------------------------------------------------------- %
global delta g gamma s sg lnpd debug rhow sig_w sig
s1=strans(s,w);
inside = delta * exp(g*(1-gamma))*exp(1/2 * (1-rhow^2) * sig_w^2 ) * exp(-gamma * (s1 - s )) .*  (1+exp(interp(s1,sg,lnpd))) .* exp((rhow * sig_w/sig - gamma)*w);
% debug(:,2)=inside';

end