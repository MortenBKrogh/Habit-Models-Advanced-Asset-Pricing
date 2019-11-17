%% Table 1 - BusinessCycles
name = string(['../Tables/BC.tex']);
if isfile(name)
delete(name);
end

%%
diary(name);
diary on

disp(['\begin{table}[H] ']);
disp(['\centering']);
disp(['\caption{Business Cycle, Simulated and historic}']);
disp(['\label{tab:BC}']);
disp(['\begin{tabular}{@{}llll@{}}']);
disp(['\toprule']);
disp(['              & \textit{Historic} & \multicolumn{2}{c}{\textit{Simulated}} \\ \midrule']);
disp(['                       &                   & $\Bar{S}$       & $\Bar{S}_{rec}$      \\ \cmidrule(l){3-4} ']);
disp(['\textit{Parameter}     &                   & ',num2str(exp(s_bar),2),'&',   num2str(exp(Rec_s_bar),2)     ,'              \\']);
disp(['Recession, \% &', num2str(rec_emp_percentage*100,4),'&',   num2str(Model_Rec*100,4),' &',   num2str(Match_Rec*100,4),'\\ \bottomrule']);
disp(['\end{tabular}']);
disp(['\end{table}']);
diary off