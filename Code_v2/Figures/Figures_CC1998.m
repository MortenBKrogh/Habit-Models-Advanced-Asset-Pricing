global Regressions
%% Figures
clear
% Figure 7 in CC1998
load('PC_Claim_workspace')

Sample = 1:500;

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
plot(S,exp(PC_ratio)/tsc,'LineWidth',1.5);% Annulized P/C-curve
hold on;
plot(S,exp(PD_ratio)/tsc,'LineWidth',1.5); % Annulized P/D-curve
ylabel('$P/C,\qquad P/D$','Interpreter','latex');
xlabel('Surplus Consumption ratio, $S$','Interpreter','latex');
xline(exp(Rec_s_bar),'--','$\bar{S}_{REC}$','Interpreter','latex');
xline(S_max,'--','$\bar{S}_{MAX}$','Interpreter','latex');
xline(0.02,'--','$\bar{S}_{2,REC}$','Interpreter','latex');
legend('PC-Ratio', 'PD-Ratio','Location','best')
hold off;
saveas(gcf,string(['../Figures/PC_PD_Ratio']),'eps2c');
%%
load('PD_Claim_workspace','elnr_pf','lnrf_pf');lnrPD = elnr_pf;
load('PC_Claim_workspace','elnr_pf');lnrPC = elnr_pf;
figure;
plot(S,lnrPC *tsc*100,'LineWidth',1.5);
hold on
plot(S,lnrPD*tsc*100,'LineWidth',1.5);
yline(mean(lnrf_pf)*tsc*100,':');
xline(exp(Rec_s_bar),'--','$\bar{S}_{REC}$','Interpreter','latex');
xline(S_max,'--','$\bar{S}_{MAX}$','Interpreter','latex');
xline(0.02,'--','$\bar{S}_{2,REC}$','Interpreter','latex');
ylabel('Expected Returns, annual percentage, $E_t ( r_{t+1} )$','Interpreter','latex');
xlabel('Surplus Consumption ratio, $S$','Interpreter','latex');
legend('Expected Return, Consumption Claim','Expected Return, Dividend Claim','Risk Free Rate','Interpreter','latex','Location','best');
saveas(gcf,string(['../Figures/ErPCPD']),'eps2c');