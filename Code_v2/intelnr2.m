function [out] = intelnr2(v)

% Fornece o valor do retorno esperado de acordo com uma distribuição normal
% ----------------------------------------------------------------------- %

global sig
out = (1/(sig*(2*pi)^(.5)))*exp(-.5*(v/sig).^2).*log(erinsd(v)).^2;
end