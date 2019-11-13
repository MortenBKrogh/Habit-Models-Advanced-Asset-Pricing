%% Figures
% Figure 7 in CC1998
figure;
subplot(2,1,1)
scatter(lnrtsim*1e2,lndctsim*1e2);title("Monthly Returns vs. consumption growth");
subplot(2,1,2)
scatter(alnrtsim_pf*1e2,alndctsim_pf*1e2);title("Annual Returns vs. consumption growth");

%% Stationary Density 
warning('off','all'); % fplot doest like the integral functions
figure;
fplot(@q_s, [min(log(S)+3) s_max+0.15]);title('Stationary Distribution of s')
hold on;
xline(s_bar, '--','$\bar{s}$', 'interpreter','latex', 'fontsize', 16);
hold on;
xline(s_max, '--','$s_{max}$', 'interpreter','latex', 'fontsize', 16);
warning('on','all');
%%
figure;
plot(S,PC_ratio/tsc,'red');title("PC vs. PD"); % Annulized P/C-curve
hold on;
plot(S,PD_ratio/tsc,'blue'); % Annulized P/D-curve
legend('PC-Ratio', 'PD-Ratio')
hold off;
%% Regression
% plot fitted vs actual
figure;
plot(y)
hold on;
plot(reg.yhat)
