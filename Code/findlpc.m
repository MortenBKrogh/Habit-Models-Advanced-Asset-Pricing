<<<<<<< HEAD
<<<<<<< HEAD
function [lnpca, ctrindx]=findlpc(sig,g,s_bar)
% This is the procedure that will calculate the fixed point P / C. %
% ----------------------------------------------------------------------- %
global s sg lnpc delta gamma phi B debug
%% S_bar
%Now we need to find the index of the value of s_bar to be used in
% graphs and other statistics
if max(sg == s_bar) == 1
ctrindx = find(sg == s_bar);
else
disp ('ERROR: The stationary value of log (S) is not in the grid');
end
%% Function value vectors P / C
lnpca = zeros(size(sg,1),1); % We are starting with P/C = 1
lnpc = lnpca;
newlnpc = lnpc;
%% Loop: find ln (P / C) from the grid of s
iter = 1;
erro = 1;
while iter < 10000 && erro > 1e-6
for i=1:size(sg,1)
s = sg(i);
if exp(-log(delta)+gamma*g-(gamma*(1-phi)-B)/2-B*(s-s_bar)) < g
disp('\t ATENÇÃO: Rf < g \n');
fprintf('value of st: %g',s);
end
% Generate the log of the variable interest rate in time
%newlnpc(i) = log(intgl1(@pdint,abs(sig)*(-8),abs(sig)*8));
newlnpc(i)=log(GaussLegendre(@pdint,abs(sig)*(-8),abs(sig)*8,40)  );
debug;
end
tv = max(abs((exp(newlnpc)-exp(lnpc))./exp(newlnpc)));
lnpc = newlnpc;
erro = max(tv);
iter = iter + 1;
end
lnpca = lnpc;
=======
=======
>>>>>>> 776d6f25bdaddaf60dfe43ccc5d8949fe61b9992
function [lnpca ctrindx]=findlpc(sig,g,s_bar)

% Este é o procedimento que irá calcular o ponto fixo P/C. %
% ----------------------------------------------------------------------- %
global s sg lnpc delta gamma phi B debug
%% S_bar
% Agora precisamos achar o index do valor de s_bar para ser utilizado em % gráficos e outras estatísticas
if max(sg == s_bar) == 1
    ctrindx = find(sg == s_bar);
else
    disp('ERRO: O valor estacionário de log(S) não está no grid');
end
%% Vetores da função valor P/C
lnpca = zeros(size(sg,1),1); % Estamos começando com P/C = 1
lnpc = lnpca;
newlnpc = lnpc;
%% Loop: achar ln(P/C) a partir do grid de s
iter = 1;
erro = 1;
while iter < 10000 && erro > 1e-6
    for i=1:size(sg,1)
        s = sg(i);
        if exp(-log(delta)+gamma*g-(gamma*(1-phi)-B)/2-B*(s-s_bar)) < g
            disp('\t ATENÇÃO: Rf < g \n');
            fprintf('valor de st: %g',s);
        end
        % Gerar o Log da taxa de juros variável no tempo
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