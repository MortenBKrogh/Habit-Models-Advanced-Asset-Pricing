function out = intercb(v)
% Integrando dos retornos esperados dos títulos públicos                  % 
% ou, da estrutura a termo.                                               % 

global sig

out = pdf('norm',v,0,sig).*log(ercbin(v));

end