function [out]=intemrs(v)

% It will return the marginal rate of substitution in such a way that it
% can be used for fixed point integration.                                %
% ----------------------------------------------------------------------- %

global sig
out = (1/(sig*(2*pi)^(.5)))*exp(-.5*(v/sig).^2).*mrsinsd(v);
end