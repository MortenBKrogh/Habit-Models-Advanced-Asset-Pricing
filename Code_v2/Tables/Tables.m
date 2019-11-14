

%% Table 1 - Calibrated Model
delete Table_1.tex
diary("Table_1.tex")
diary on

disp('\begin{table}[H]')
disp('\centering')
disp('\begin{threeparttable}[b]')
disp('\caption{Parameters of the model}')
disp('\label{tab:ModelCalib}')
disp('\begin{tabular}{@{}ll@{\hspace{1.5cm}}ll@{}}')
disp('\toprule')
disp(' & Parameter                              & Notation         & Value    \\ \midrule') 
disp('\multicolumn{4}{l}{\textit{Calibrated}}                                 \\')
disp([' & Mean consumption growth                & $g$            & $', num2str(Edc_pf), '$ \\'])
disp([' & Standard deviation of $\Delta c_t$     & $\sigma$         & $', num2str(Stdc_pf),'$ \\'])
disp([' & Standard deviation of $\Delta d_t$     & $\sigma_w$       & $', num2str(Stdc_pf), '$ \\'])
disp([' & Log risk-free rate                     & $r^f$            & $', num2str(Erfinterp_pf),'$ \\'])
disp([' & Persistence parameter                  & $\phi$           & $', num2str(phi),'$ \\'])
disp([' \multicolumn{4}{l}{\textit{Assumed}}                                   \\'])
disp([' & Coefficient of Risk Aversion           & $\gamma$         & $', num2str(gamma),'$ \\'])
disp([' & Correlation dividends/consumption      & $\rho$           & $', num2str(phi),'$ \\'])
disp(['\multicolumn{4}{l}{\textit{Implied}}                                    \\'])
disp([' & Subjective discount factor             & $\delta$         & $', num2str(delta),'$ \\'])
disp(['& Steady-state surplus consumption ratio & $\Bar{S}$        & $', num2str(S_bar),'$ \\'])
disp([' & Maximum surplus consumption ratio      & $S_{\text{max}}$ & $', num2str(S_max),'$ \\ \bottomrule'])
disp('\end{tabular}')
disp('\begin{tablenotes}')
disp('\footnotesize{\item [1] All relevant parameters are annualized')
disp('              \item [2] Calibrated parameters are estimated from data, assumed are chosen arbitrarily on the grounds of existing literature, while implied parameters are calculated from the calibrated/assumed parameters.}')
disp('\end{tablenotes}')
disp('\end{threeparttable}')
disp('\end{table}')

diary off

%% Table 2 - Data Properties
delete Table_2.tex
diary("Table_2.tex")
diary on
disp(['\begin{table}[H]'])
disp(['\centering'])
disp(['\caption{Data Properties}'])
disp(['\label{tab:Data_props}'])
disp(['\begin{tabular}{@{}l@{\hspace{1.5cm}}l@{\hspace{1.5cm}}l@{}}'])
disp(['\toprule'])
disp([' & \textit{Simulated} & \textit{Historic} \\ \midrule'])
disp(['$\mathbb{E}\left[r_t- r^f_t\right]$& $', num2str(Eexrettinterp_pf),'$           & $', num2str(0.0927),'$          \\'])
disp(['$\sigma\left(r_t - r^f_t  \right)$ & $', num2str(Stdexrettinterp_pf),'$           & $', num2str(0.1670),'$          \\'])
disp(['$\mathbb{E}\left[r_t- r^f_t\right] / \sigma\left(r_t - r^f_t,\right)$ & $', num2str(Shprinterp_pf),'$ & $', num2str(0.5548),'$  \\ \bottomrule'])
disp(['\end{tabular}'])
disp(['\end{table}'])
diary off


%% Table 3 - Simulated Moments
delete Table_3.tex
diary("Table_3.tex")
diary on
disp(['\begin{table}[H]'])
disp(['\centering'])
disp(['\caption{Simulated Moments}'])
disp(['\label{tab:MMoomme}'])
disp(['\begin{tabular}{@{}llllllllll@{}}'])
disp(['\toprule '])
disp([' & $\mathbb{E}\Delta d$ & $\sigma_{\Delta d}$ & $\mathbb{E}r^f$ & $\mathbb{E}r^m/\sigma _{r^m}$ & $\mathbb{E}R^m/\sigma _{R^m}$ & $\mathbb{E}r^m$ & $\sigma_{r^m}$ & $\mathbb{E}d-p$ & $\sigma_{d-p}$  \\ '])
disp(['\midrule '])
disp(['\multicolumn{10}{l}{$P/D$}\\'])
disp([' & ', num2str(0.0117),' & ', num2str(0.1026),' & ', num2str(0.0109),' & ', num2str(0.2424),' & ', num2str(0.3147),' & ', num2str(0.0402),' & ', num2str(0.1657),' & ', num2str(3.2545),' & ', num2str(0.2070),' \\ '])
disp(['\multicolumn{10}{l}{$P/C$}\\'])
disp(['& ', num2str(0.0135),' & ', num2str(0.0124),' & ', num2str(0.0109),' & ', num2str(0.3853),' & ', num2str(0.4204),' & ', num2str(0.0372),'& ', num2str(0.0966),' & ', num2str(3.3797),' & ', num2str(0.1864),' \\ '])
disp(['\bottomrule '])
disp(['\end{tabular}'])
disp(['\end{table}'])
diary off

 