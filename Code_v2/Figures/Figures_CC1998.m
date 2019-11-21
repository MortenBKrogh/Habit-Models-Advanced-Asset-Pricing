global Regressions
%% Figures
clear
% Figure 7 in CC1998
load('PD_Claim_workspace')

Sample = 100:1000;

figure;
subplot(2,1,1)
scatter(lnrtsim(Sample)*1e2,lndctsim(Sample)*1e2);title("Monthly Returns vs. consumption growth");
subplot(2,1,2)
scatter(alnrtsim_pf(Sample)*1e2,alndctsim_pf(Sample)*1e2);title("Annual Returns vs. consumption growth");
hold off
%saveas(gcf,string(['Figures/Figure_7_CC_1998_Calib_', num2str(calib),'_PD_', num2str(PD_Claim), '.eps']),'eps2c');


%%
load('PD_Claim_workspace','output_lnpda','S','tsc');PD_ratio = output_lnpda;
load('PC_Claim_workspace','output_lnpca');PC_ratio = output_lnpca;
%%
figure;
plot(S,PC_ratio/tsc);% Annulized P/C-curve
hold on;
plot(S,PD_ratio/tsc); % Annulized P/D-curve
ylabel('$P/C,\qquad P/D$','Interpreter','latex');
legend('PC-Ratio', 'PD-Ratio','Location','northwest')
xlabel('Surplus Consumption ratio, $S_t$','Interpreter','latex');
hold off;
saveas(gcf,string(['../Figures/PC_PD_Ratio']),'eps2c');
%%
load('PD_Claim_workspace','elnr_pf');lnrPD = elnr_pf;
load('PC_Claim_workspace','elnr_pf');lnrPC = elnr_pf;
figure;
plot(S,lnrPC*tsc*100);
hold on
plot(S,lnrPD*tsc*100);
xline(exp(Rec_s_bar),'--','$\bar{S}_{REC}$','Interpreter','latex');
xline(S_max,'--','$\bar{S}_{MAX}$','Interpreter','latex');
ylabel('Expected Returns, annual percentage, $E_t ( r_{t+1} )$','Interpreter','latex');
xlabel('Surplus Consumption ratio, $S_t$','Interpreter','latex');
legend('Expected Return, Consumption Claim','Expected Return, Dividend Claim','Interpreter','latex');
saveas(gcf,string(['../Figures/ErPCPD']),'eps2c');

%% Regression
% plot fitted vs actual
if Regressions == 1
figure;
plot(y)
hold on
plot(reg.yhat)
saveas(gcf,string(['Figures/Reg_Fitted_vs_Actual', num2str(calib),'_PD_', num2str(PD_Claim), '.eps']),'eps2c');
end