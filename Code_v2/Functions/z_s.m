function [z] = z_s(s)
% Functional Z of s, for the theoretical density of S
global gamma S_bar sig
lnZ =  -gamma .* S_bar^2 .* ((S_bar^(-2)-1)./(lambda(s))+3.*lambda(s)+(lambda(s).^2)./2) - (gamma .* (3 .* S_bar^2-1) +2) .* log(lambda(s))-2*log(sig);
z = exp(lnZ);
end

