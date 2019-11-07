function [inside] = pdivint(w)
% Will create a new normal density according to the innovations of w {t + 1}
% to integrate P / D functional.
% ----------------------------------------------------------------------- % 
global sig_w debug
inside = (1/(sig_w*(2*pi)^(.5)))*exp(-.5*(w/sig_w).^2).*pdivsmotor(w);
debug(:,3)=inside'; 
end
