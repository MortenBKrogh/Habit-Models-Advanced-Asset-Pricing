%% Figures
% Figure 7 in CC1998
figure;
subplot(2,1,1)
scatter(lnrtsim*1e2,lndctsim*1e2);title("Monthly Returns vs. consumption growth");
subplot(2,1,2)
scatter(alnrtsim_pf*1e2,alndctsim_pf*1e2);title("Annual Returns vs. consumption growth");

%% Stationary Density 
figure;
fplot(@q_s, [min(log(S)) s_max]);title('Stationary Distribution of s')