function [er]=erinsd(v)
% Procedure for calculating consumption claim returns.                    % 
% ----------------------------------------------------------------------- % 

global sg s lnpca g 
er=((1+exp(interp(strans(s,v),sg,lnpca)))'./(ones(size(v))...
*exp(interp(s,sg,lnpca)))).*exp(g+v); 

end