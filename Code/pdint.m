function [inside] = pdint(v)
% Irá criar uma nova densidade normal de acordo com as inovações de v{t+1} 
% para integrar o funcional de P/C.
% ----------------------------------------------------------------------- % 
global sig debug
inside = (1/(sig*(2*pi)^(.5)))*exp(-.5*(v/sig).^2).*pdmotor(v);
debug(:,3)=inside'; 
end
