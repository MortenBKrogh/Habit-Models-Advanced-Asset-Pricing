function [inside] = pdivint(w)
% Will create a new normal density according to the innovations of w {t + 1}
% to integrate P / D functional.
% ----------------------------------------------------------------------- % 
global sig debug
%inside = 1/(sig* (2*pi)^(1/2) ) .* exp(-(1/2) .* (w/sig).^2) .*pdivsmotor(w);
inside =(1/(sig* (2*pi)^(.5 ) ))*exp(-.5*(w/sig).^2).*pdivsmotor(w);
% debug(:,3)=inside'; 
end
