A = struct2cell(regs);
name = string(['Tables/RegressionTable.tex']);
if isfile(name)
delete(name);
end
%%
diary(name);
diary on

disp(['\begin{table}[H]']);
disp(['\centering   ']);
disp(['  \caption{Regressions}           ']);
disp(['  \label{tab:regress}     ']);
disp(['  \begin{threeparttable}']);    
disp(['\begin{tabular}{@{\hspace{5pt}}l@{\hspace{5pt}}cccc} ']);
disp(['\toprule ']);
disp([' & \multicolumn{4}{c}{\textit{Dependent variable:}} \\ ']);
disp([' & \multicolumn{4}{c}{$\left(r_{t+1}-r^f\right)$} \\ ']);
disp([' \cmidrule(rr){2-5}']);
disp([' & (1)   &   (2) &   (3) &  (4)\\ ']);
disp(['\midrule  ']);
disp(['\\[-2.1ex] $\left( p_t - c_t \right)_{REC}$ & $-$',num2str(abs(A{5,1,3}(2)),4),'& &  \\ ']);
disp(['  & (',num2str(abs(A{10,1,3}(2)),2),') & & & \\ ']);
disp([' \addlinespace ']);
disp(['  $\left( p_t - c_t \right)_{EXP}$ & $-$', num2str(abs(A{5,1,3}(3)),4),' & &  \\ ']);
disp(['  & (', num2str(abs(A{10,1,3}(3)),2),') & & &\\ ']);
disp([' \addlinespace ']);
disp([' $p_t - c_t$ &  & $-$',num2str(abs(A{5,1,4}(2)),4),' & & \\']);
disp([' & & (',num2str(abs(A{10,1,4}(2)),2),') \\']);
disp([' \addlinespace ']);
disp(['  $\left( p_t - d_t \right)_{REC}$ & & & $-$',num2str(abs(A{5,1,1}(2)),4),'&  \\ ']);
disp(['  & & & (',num2str(abs(A{10,1,1}(2)),2),')    &\\ ']);
disp([' \addlinespace ']);
disp(['  $\left( p_t - d_t \right)_{EXP}$ & & & $-$', num2str(abs(A{5,1,1}(3)),4),'&  \\ ']);
disp(['  &  & & (',num2str(abs(A{10,1,1}(3)),2), ') &\\ ']);
disp([' \addlinespace ']);
disp([' $p_t - d_t$ & & & & $-$',num2str(abs(A{5,1,2}(2)),4),' \\']);
disp([' & & & &  (',num2str(abs(A{10,1,2}(2)),2),')  \\']);
disp([' \addlinespace ']);
disp([' Constant &',num2str(abs(A{5,1,3}(1)),4),' &' num2str(abs(A{5,1,4}(1)),4),' &', num2str(abs(A{5,1,1}(1)),4),' &',num2str(abs(A{5,1,2}(1)),4),' \\ ']);
disp(['  &(',num2str(abs(A{10,1,3}(1)),2),') &(', num2str(abs(A{10,1,4}(1)),2),') &(', num2str(abs(A{10,1,1}(1)),2),') &(',num2str(abs(A{10,1,2}(1)),2),') \\ ']);
disp([' \addlinespace ']);
disp(['\midrule  ']);
disp(['Observations & ',num2str(abs(A{3,1,3})),' & ',num2str(abs(A{3,1,4})),' & ',num2str(abs(A{3,1,1})),' &',num2str(abs(A{3,1,2})),'\\ ']);
disp(['R$^{2}$ &',num2str(abs(A{11,1,3}),2),' & ',num2str(abs(A{11,1,4}),2),' & ',num2str(abs(A{11,1,1}),2),' &',num2str(abs(A{11,1,2}),2),'\\ ']);
disp(['Residual Std. Error &',num2str(abs(A{8,1,3}),2),' & ',num2str(abs(A{8,1,4}),2),' &',num2str(abs(A{8,1,1}),2),' & ',num2str(abs(A{8,1,2}),2),'  \\ ']);
disp(['\bottomrule ']);
disp(['\end{tabular} ']);
disp(['\begin{tablenotes}']);
disp(['\footnotesize{']);
disp(['\item[1] Brackets below estimates contains Newey-West corrected standard errors. ']);
disp(['\item[2] Regressions on ',num2str(abs(A{3,1,3})),' years of simulated data.']);
disp(['\item[3] EXP (REC) denotes expansion (recession)']);
disp(['}']);
disp(['\end{tablenotes}']);
disp(['\end{threeparttable}']);
disp(['\end{table} ']);
diary off
