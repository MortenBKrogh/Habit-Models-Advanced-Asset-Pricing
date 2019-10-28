function [inside] = pdmotor(v)
% Procedimento para ser usado na integração numérica quando for se calcular 
% o ponto fixo. Tem como argumento apenas v ~ N(0,sig). Retorna o VALOR de
% P/C para cada iteração sobre o valor corrente de s{t} em cada iteração.
% ----------------------------------------------------------------------- %
global delta g gamma s sg lnpc debug
s1=strans(s,v);
inside = delta*exp(g*(1-gamma))*exp(-gamma*(s1-s)).*...
(1+exp(interp(s1,sg,lnpc)))'.*exp((1-gamma)*v);
debug(:,2)=inside';
end