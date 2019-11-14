function [out] = inter2(v)

% Expected variance returns consumption claim
% ----------------------------------------------------------------------- %

global sig
out = (1/(sig*(2*pi)^(.5)))*exp(-.5*(v/sig).^2).*erinsd(v).^2;
end
