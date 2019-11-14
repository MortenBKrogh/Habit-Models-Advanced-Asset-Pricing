function [inside] = internorm(v)
global sig
% ------------------------------------------------------------------------%
% Density of a normal                                                     %
% ------------------------------------------------------------------------%
inside =  1/((2*pi)^(1/2) .* sig) .* exp(-(v .^2)/(2 * sig ^2));
end