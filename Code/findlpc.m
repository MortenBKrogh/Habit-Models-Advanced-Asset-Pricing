function [lnpca ctrindx]=findlpc(sig,g,s_bar)

% Este � o procedimento que ir� calcular o ponto fixo P/C.                %
% ----------------------------------------------------------------------- %

global s sg lnpc delta gamma phi B debug

%% S_bar
% Agora precisamos achar o index do valor de s_bar para ser utilizado em
% gr�ficos e outras estat�sticas

if max(sg == s_bar) == 1

    ctrindx = find(sg == s_bar);
else
    disp('ERRO: O valor estacion�rio de log(S) n�o est� no grid');
end

%% Vetores da fun��o valor P/C
lnpca = zeros(size(sg,1),1); % Estamos come�ando com P/C = 1 lnpc = lnpca;
newlnpc = lnpc;

%% Loop: achar ln(P/C) a partir do grid de s iter = 1;

erro = 1;

while iter < 10000 && erro > 1e-6
    for i=1:size(sg,1) s = sg(i);
        if exp(-log(delta)+gamma*g-(gamma*(1-phi)-B)/2-B*(s-s_bar)) < g disp('\t ATEN��O: Rf < g \n');
            fprintf('valor de st: %g',s);
        end
        % Gerar o Log da taxa de juros vari�vel no tempo
        newlnpc(i) = log(intgl1(@pdint,abs(sig)*(-8),abs(sig)*8));
        debug;
    end
    tv = max(abs((exp(newlnpc)-exp(lnpc))./exp(newlnpc)));
    lnpc = newlnpc;
    erro = max(tv);
    iter = iter + 1;
end
lnpca = lnpc;
=======
function [lnpca ctrindx]=findlpc(sig,g,s_bar)

% Este � o procedimento que ir� calcular o ponto fixo P/C. %
% ----------------------------------------------------------------------- %
global s sg lnpc delta gamma phi B debug
%% S_bar
% Agora precisamos achar o index do valor de s_bar para ser utilizado em % gr�ficos e outras estat�sticas
if max(sg == s_bar) == 1
    ctrindx = find(sg == s_bar);
else
    disp('ERRO: O valor estacion�rio de log(S) n�o est� no grid');
end
%% Vetores da fun��o valor P/C
lnpca = zeros(size(sg,1),1); % Estamos come�ando com P/C = 1
lnpc = lnpca;
newlnpc = lnpc;
%% Loop: achar ln(P/C) a partir do grid de s
iter = 1;
erro = 1;
while iter < 10000 && erro > 1e-6
    for i=1:size(sg,1)
        s = sg(i);
        if exp(-log(delta)+gamma*g-(gamma*(1-phi)-B)/2-B*(s-s_bar)) < g
            disp('\t ATEN��O: Rf < g \n');
            fprintf('valor de st: %g',s);
        end
        % Gerar o Log da taxa de juros vari�vel no tempo
        newlnpc(i) = log(intgl1(@pdint,abs(sig)*(-8),abs(sig)*8));
        debug;
    end
    tv = max(abs((exp(newlnpc)-exp(lnpc))./exp(newlnpc)));
    lnpc = newlnpc;
    erro = max(tv);
    iter = iter + 1;
end
lnpca = lnpc;
end
