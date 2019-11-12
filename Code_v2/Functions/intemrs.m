function [out]=intemrs(v)

% Retornará a taxa marginal de substituição de tal forma que possa ser    %
% utilizada na integração do ponto fixo.                                  %
% ----------------------------------------------------------------------- %

global sig
out = (1/(sig*(2*pi)^(.5)))*exp(-.5*(v/sig).^2).*mrsinsd(v);
end