
name = string(['../Tables/SimulatedMom.tex']);
if isfile(name)
delete(name);
end
%%
disp(['\begin{table}[H]']);
disp(['\centering']);
disp(['\caption{Simulated Moments}']);
disp(['\label{tab:simmom}']);
disp(['\begin{tabular}{@{}lllll@{}}']);
disp(['\toprule']);
disp(['Statistic                                               & \makecell{Consumption \\ Claim} & \makecell{Dividend \\ Claim} & \makecell{CC99-Calibration \\ Consumption Claim} & \makecell{CC99-Calibration\\ Dividend Claim} \\ \midrule']);
disp(['$\mathbb{E}\left(\Delta c \right)$                      &',num2str(PC_edc),'&',num2str(PD_edc),'&',num2str(CC_PC_edc),'&',num2str(CC_PD_edc),'\\']);
disp(['$\sigma\left(\Delta c \right)$                          &',num2str(PC_stdc),'&',num2str(PD_stdc),'&',num2str(CC_PC_stdc),'&',num2str(CC_PD_stdc),'\\']);
disp(['$\mathbb{E}r^f$                                         &',num2str(PC_rfr),'&',num2str(PD_rfr),'&',num2str(CC_PC_rfr),'&',num2str(CC_PD_rfr),'\\']);
disp(['$\mathbb{E}\left(r-f^f\right)/\sigma\left(r-r^f\right)$ &',num2str(PC_shrp),'&',num2str(PD_shrp),'& ',num2str(CC_PC_shrp),'                                       &         ',num2str(CC_PD_shrp),'                            \\']);
disp(['$\mathbb{E}\left(R-R^f\right)/\sigma\left(R-R^f\right)$ & ',num2str(PC_SHRP),'   & ',num2str(PD_SHRP),'               &  ',num2str(CC_PC_SHRP),'                                      &      ',num2str(CC_PD_SHRP),'                               \\']);
disp(['$\mathbb{E}\left(r-r^f\right)$                          &  ',num2str(PC_Eer),'                 &      ',num2str(PD_Eer),'          &         ',num2str(CC_PC_Eer),'                               &        ',num2str(CC_PD_Eer),'                             \\']);
disp(['$\sigma\left(r-r^f\right)$                              &     ',num2str(PC_Stder),'              &       ',num2str(PD_Stder),'           &   ',num2str(CC_PC_Stder),'                                       &                  ',num2str(PD_Stder),'                     \\']);
disp(['$\mathbb{E}\left(p-c\right)$                          &     ',num2str(PC_Epd),'                &      ',num2str(PD_Epd),'            &         ',num2str(CC_PC_Epd),'                                 &             ',num2str(CC_PD_Epd),'                          \\']);
disp(['$\sigma\left(p-c\right)$                              &       ',num2str(PC_StdPD),'            &      ',num2str(PD_StdPD),'          &   ',num2str(CC_PC_StdPD),'                                     & ',num2str(CC_PD_StdPD),'                                    \\ \bottomrule']);
disp(['\end{tabular}']);
disp(['\end{table}']);
