function [out]=mrsinsd(v)
% Retorna a taxa marginal de substitui??o intertemporal no modelo. % 
% ----------------------------------------------------------------------- 
global delta g gamma s
out = delta*exp(-gamma*g)*exp(-gamma*v).*exp(-gamma*(strans(s,v)-s));
end