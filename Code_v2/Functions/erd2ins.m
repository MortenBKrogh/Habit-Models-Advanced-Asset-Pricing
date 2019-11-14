function [out] = erd2ins(v)

% Expected variance returns consumption claim
% ----------------------------------------------------------------------- %

global sig g rhow sig_w lnpca sg
out = (1 + exp(interp(strans(s,v), sg, lnpca)))...
                            ./exp(interp(s,sg,lnpca)) .* exp(g) .*...
                exp( rhow.* sig_w./ sig.*v).*exp((1- rhow^2) * sig_w^2)^2;

end