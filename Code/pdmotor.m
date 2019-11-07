function [inside] = pdmotor(v)
% Procedure to be used for numeric integration when calculating
% the fixed point. It has as argument only v ~ N (0, sig). Returns VALUE of
% P / C for each iteration over the current value of s {t} in each iteration.
% ----------------------------------------------------------------------- %
global delta g gamma s sg lnpc debug
s1=strans(s,v);
inside = delta*exp(g*(1-gamma))*exp(-gamma*(s1-s)).*...
(1+exp(interp(s1,sg,lnpc)))'.*exp((1-gamma)*v);
debug(:,2)=inside';
end