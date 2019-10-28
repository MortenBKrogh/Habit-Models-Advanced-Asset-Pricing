function out = intpcb(v)

% Função que fornece o preço dos bonds para cada maturidade               %

global s sg lnpcb sig

out = pdf('norm',v,0,sig).*mrsinsd(v).*...
    exp(interp(strans(s,v),sg,lnpcb(:,size(lnpcb,2))))';

end