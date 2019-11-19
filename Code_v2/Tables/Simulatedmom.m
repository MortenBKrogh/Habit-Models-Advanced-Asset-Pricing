
name = string(['../Tables/SimulatedMom.tex']);
if isfile(name)
delete(name);
end
%%
diary(name);
diary on

disp(['\begin{table}[H]']);
disp(['\centering']);
disp(['\caption{Simulated Moments}']);
disp(['\label{tab:simmom}']);
disp(['\begin{tabular}{@{}lcccc@{}}']);
disp(['\toprule']);
disp(['Statistic                                               & \makecell{Consumption \\ Claim} & \makecell{Dividend \\ Claim} & \makecell{CC99-Calibration \\ Consumption Claim} & \makecell{CC99-Calibration\\ Dividend Claim} \\ \midrule']);
disp(['$\mathbb{E}\left(\Delta c \right)$                      &',num2str(PC_edc,3),'&',num2str(PD_edc,3),'&',num2str(CC_PC_edc,3),'&',num2str(CC_PD_edc,3),'\\']);
disp(['$\sigma\left(\Delta c \right)$                          &',num2str(PC_stdc,3),'&',num2str(PD_stdc,3),'&',num2str(CC_PC_stdc,3),'&',num2str(CC_PD_stdc,3),'\\']);
disp(['$\mathbb{E}r^f$                                         &',num2str(PC_rfr,3),'&',num2str(PD_rfr,3),'&',num2str(CC_PC_rfr,3),'&',num2str(CC_PD_rfr,3),'\\']);
disp(['$\mathbb{E}\left(r-f^f\right)/\sigma\left(r-r^f\right)$ &',num2str(PC_shrp,3),'&',num2str(PD_shrp,3),'& ',num2str(CC_PC_shrp,3),'                                       &         ',num2str(CC_PD_shrp,3),'                            \\']);
disp(['$\mathbb{E}\left(R-R^f\right)/\sigma\left(R-R^f\right)$ & ',num2str(PC_SHRP,3),'   & ',num2str(PD_SHRP,3),'               &  ',num2str(CC_PC_SHRP,3),'                                      &      ',num2str(CC_PD_SHRP,3),'                               \\']);
disp(['$\mathbb{E}\left(r-r^f\right)$                          &  ',num2str(PC_Eer,3),'                 &      ',num2str(PD_Eer,3),'          &         ',num2str(CC_PC_Eer,3),'                               &        ',num2str(CC_PD_Eer,3),'                             \\']);
disp(['$\sigma\left(r-r^f\right)$                              &     ',num2str(PC_Stder,3),'              &       ',num2str(PD_Stder,3),'           &   ',num2str(CC_PC_Stder,3),'                                       &                  ',num2str(PD_Stder,3),'                     \\']);
disp(['$\mathbb{E}\left(p-c\right)$                          &     ',num2str(PC_Epd,3),'                &      ',num2str(PD_Epd,3),'            &         ',num2str(CC_PC_Epd,3),'                                 &             ',num2str(CC_PD_Epd,3),'                          \\']);
disp(['$\sigma\left(p-c\right)$                              &       ',num2str(PC_StdPD,3),'            &      ',num2str(PD_StdPD,3),'          &   ',num2str(CC_PC_StdPD,3),'                                     & ',num2str(CC_PD_StdPD,3),'                                    \\ \bottomrule']);
disp(['\end{tabular}']);
disp(['\end{table}']);
diary off