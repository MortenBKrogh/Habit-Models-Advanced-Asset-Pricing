function out = ercbin(v)
% Retornos esperados da estrutura a termo %
global s sg lnpcb matur
out = exp(interp(strans(s,v),sg,lnpcb(:,matur-1))) ./ ...
    exp(interp(s,sg,lnpcb(:,matur)));
end