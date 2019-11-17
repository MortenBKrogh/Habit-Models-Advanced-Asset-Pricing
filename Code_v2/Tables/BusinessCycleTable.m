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
disp(['\begin{tabular}{@{\hspace{2mm}}ll@{\hspace{5mm}}ll@{\hspace{5mm}}}']);
disp(['\toprule']);
disp(['                       & \multicolumn{2}{c}{\textit{Simulated}} & \textit{Historic}  \\ \midrule']);
disp(['                                          & $\Bar{S}$       & $\Bar{S}_{rec}$ &     \\ \cmidrule(l){2-3} ']);
disp(['\textit{Parameter}                        & ',num2str(exp(s_bar),2),'&',   num2str(exp(Rec_s_bar),2)     ,'&              \\']);
disp(['\textit{Recession}, \% &',   num2str(Model_Rec*100,4),' &',   num2str(Match_Rec*100,4),'&', num2str(rec_emp_percentage*100,4),'\\ \bottomrule']);
disp(['\end{tabular}']);
disp(['\end{table}']);
diary off