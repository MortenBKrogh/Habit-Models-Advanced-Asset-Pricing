function [inside] = erdinsd(v)
global g s sg rhow sig sig_w
% ------------------------------------------------------------------------%
% Numerical integration expected returns P/D                              %
% ------------------------------------------------------------------------%
inside = (1 + exp(interp(strans(s,v),sg,lnpd)))...
                            ./exp(interp( s, sg, lnpd1)) .* exp(g) .*...
                exp(rhow.* sig_w./ sig .* v).* exp(1/2 * (1- rhow^2)* sig_w^2);
end