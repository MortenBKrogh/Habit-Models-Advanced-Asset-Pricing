function [er elnr sdr sdlnr lnrf lnrf1 lny elnrcb sdlnrcb slpmv] = finders(sg)

% Procedimento que calcula os retornos esperados dos ativos de consumo. ele
% nos fornecer? E(R), SD(R), lnrf, a curvatura da fronteira m?dia vari?ncia
% dada pela vari?vel (slpmv) e integrar? a raz?o de Sharpe do vetor de    %
% ativos de consumo.                                                      %
% ----------------------------------------------------------------------- %

global g gamma sig phi s maxcb s_bar delta tsc lnpcb matur

% Inclina??o da Fronteira de M?dia-Vari?ncia                              %
%                                                                         %
% slpmv = (exp((gamma*sig)^2.*(1+lambda(sg)).^2)-1).^.5                   %
% ----------------------------------------------------------------------- %

slpmv = (exp((gamma*sig)^2*(1+lambda(sg)).^2)-1).^(0.5);

%% Estrutura termo da taxa de juros dada por:                             %
%                                                                         %
% lnrf = -ln(delta) + gamma*g -                                           %
% gamma*(1-phi)*(s{t}-s_bar)-.5((gamma*sig)^2)*(1+lambda(s{t}))^2
% - B*(sg - s_bar) --> Vari?vel no estado
% ----------------------------------------------------------------------- %

lnrf = -log(delta) + gamma*g - gamma*(1-phi)*(sg-s_bar)...
    - 0.5*(gamma*sig*(1+lambda(sg))).^2;

%% Bonds
% Matriz de todos os pre?os dos t?tulos. Sua dimens?o ser?                
% N(sg) x (maxcb*tsc)
lnpcb = [];
lnpcb(:,1) = -lnrf;

lnp = zeros(size(sg,1),1);

for j = 2:(maxcb*tsc)
    for i = 1:length(sg)
        s = sg(i);
        lnp(i) = log(GaussLegendre(@intpcb,abs(sig)*(-8),abs(sig)*8,40));
        %      lnp(i) = log(GaussLegendre(@intpcb,abs(sig)*(-8),abs(sig)*8));
        % newlnpc(i) = log(GaussLegendre(     @pdint,abs(sig)*(-8),abs(sig)*8));
        % newlnpc(i)=log(GaussLegendre(@pdint,abs(sig)*(-8),abs(sig)*8,40)  );
    end
    lnpcb = cat(2,lnpcb,lnp);
end

% Yields
lny = - ...
    lnpcb./kron(ones(size(sg,1),1),linspace(1/tsc,(maxcb*tsc)/tsc,(maxcb*tsc)));

%% Retornos e Desvios Padr?o Esperados                                    %
% ----------------------------------------------------------------------- %

lnrf1 = zeros(size(sg,1),1);

er = zeros(size(sg,1),1);
elnr = zeros(size(sg,1),1);
sdr = zeros(size(sg,1),1);
sdlnr = zeros(size(sg,1),1);
elnrcb = zeros(size(sg,1),maxcb*tsc);       % zero-coupon bonds %
sdlnrcb = zeros(size(sg,1),maxcb*tsc);

for i=1:size(sg,1)
    s = sg(i);
    
    
    lnrf1(i)= - log(GaussLegendre(@intemrs,abs(sig)*(-8),abs(sig)*8,40));
    er(i)= GaussLegendre(@inter,abs(sig)*(-8),abs(sig)*8,40);
    elnr(i)= GaussLegendre(@intelnr,abs(sig)*(-8),abs(sig)*8,40);
    sdr(i) = GaussLegendre(@inter2,abs(sig)*(-8),abs(sig)*8,40);
    sdr(i) = (sdr(i) - er(i).^2).^(.5);
    sdlnr(i) = GaussLegendre(@intelnr2,abs(sig)*(-8),abs(sig)*8,40);
    
    % Bonds
    matur = maxcb*tsc; elnrcb(i,1) = lnrf(i);
    
    for k = 2:matur
        elnrcb(i,k) = GaussLegendre(@intelnrcb,abs(sig)*(-8),abs(sig)*8, 40); 
        sdlnrcb(i,k) = GaussLegendre(@intelnr2,abs(sig)*(-8),abs(sig)*8, 40); 
        sdlnrcb(i,k) = (sdlnrcb(i,k) - elnrcb(i,k).^2).^(.5);
    end
    
end
end


