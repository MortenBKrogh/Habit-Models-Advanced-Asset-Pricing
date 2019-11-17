% Monthly
h = 1;
load('PD_Claim_workspace','lnrtsim','lnpctsim','stsim');
stsim = stsim(2:end);
lnpctsim = lnpctsim(2:end);
rec_sim_ss = NaN(length(stsim), 1);
for i = 1:length(stsim)
    if stsim(i) < Rec_s_bar
        rec_sim_ss(i) = 1;
    else
        rec_sim_ss(i) = 0;
    end
end
rec_sim_ss_percentage = sum(rec_sim_ss(:)==1) / length(rec_sim_ss);


PC_regress_m = lnpctsim;                % PD/PC
rfr  = Erfinterp_pf/12;                 % Risk free rate
rets = lnrtsim - rfr;                   % Excess Returns

y   = rets(1+h:end,1);
x   = [ones(length(rets(1:end-h,:)), 1),  ...              % vector of Ones
    rec_sim_ss(1:end-h,:) .* PC_regress_m(1:end-h,1), ...  %    I_rec_t *PD_t
    (1-rec_sim_ss(1:end-h,:)) .* PC_regress_m(1:end-h,1)]; % (1-I_rec_t)*PD_t
regPCrec = nwest(y,x,0);