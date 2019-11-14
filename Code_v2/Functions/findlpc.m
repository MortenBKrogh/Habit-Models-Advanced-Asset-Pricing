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
%% Function value vectors P / C or P / D
lnpca = zeros(size(sg,1),1); % We are starting with PC or PD = 1
lnpc = lnpca;
newlnpc = lnpc;
%% Loop: find ln (P / C) from the grid of s
iter = 1;
erro = 1;
while iter < 10000 && erro > 1e-6
for i=1:size(sg,1)
s = sg(i);
if exp(-log(delta)+gamma*g-(gamma*(1-phi)-B)/2-B*(s-s_bar)) < g
disp('\t Attention: Rf < g \n');
fprintf('value of st: %g',s);
end
% Generate the log of the variable interest rate in time
newlnpc(i)=log(GaussLegendre(@pdint,abs(sig)*(-8),abs(sig)*8,40)  );
debug;
end
tv = max(abs((exp(newlnpc)-exp(lnpc))./exp(newlnpc)));
lnpc = newlnpc;
erro = max(tv);
iter = iter + 1;
end
lnpca = lnpc;
end