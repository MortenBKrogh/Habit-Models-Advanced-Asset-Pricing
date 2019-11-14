function [inside] = pdint(v)
% Will create a new normal density according to the innovations of v {t + 1}
% to integrate P / C functional.
% ----------------------------------------------------------------------- % 
global sig debug
inside = (1/(sig*(2*pi)^(.5)))*exp(-.5*(v/sig).^2).*pdmotor(v);
debug(:,3)=inside'; 
end
