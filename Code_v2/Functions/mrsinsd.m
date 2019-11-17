function [out]=mrsinsd(v)
% Returns the marginal rate of intertemporal substitution in the model. % 
% ----------------------------------------------------------------------- 
global delta g gamma s
out = delta*exp(-gamma*g)*exp(-gamma*v).*exp(-gamma*(strans(s,v)-s));
end