load('PD_Claim_workspace','s_bar','s_max',...
        'verd','S_bar','sig','gamma','S','stsim','lnrtsim','lnpctsim','Erfinterp_pf');
    Erfinterp_pf = Erfinterp_pf./12;
    PD_regress   = lnpctsim(2:end,1);             % PD
    load('PC_Claim_workspace','lnrtsim','lnpctsim')
    alnrtsim_pf = lnrtsim;
    PC_regress  = lnpctsim(2:end,1);              % PC
    h=1;
rfr  = Erfinterp_pf;                    % Risk free rate
rets = alnrtsim_pf - rfr;               % Excess Returns
h    =  1;                              % Forecast Horizon 0 = in-sample regression
y   = rets(1+h:end,1);                  % Regressand 
%% Regressions 2
rec_sim_02 = zeros(size(stsim,1),1);
for i = 1:size(stsim,1)
    if stsim(i) < log(0.02)
        rec_sim_02(i) = 1;
    else
        rec_sim_02(i) = 0;
    end 
end
if annual == 0
rec_sim_02 = rec_sim_02(2:end);
end
%%
lower_Sbar = 1;
s_bar_2 = log(0.02);
x   = [ones(length(rets(1:end-h,:)), 1),  ...            
    rec_sim_02(1:end-h,:) .* PD_regress(1:end-h,1), ...  
    (1-rec_sim_02(1:end-h,:)) .* PD_regress(1:end-h,1)]; 
regPDrec = nwest(y,x,0); %% Full BC <- PD

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
    rec_sim_02(1:end-h,:) .* PC_regress(1:end-h,1), ...  
    (1-rec_sim_02(1:end-h,:)) .* PC_regress(1:end-h,1)]; 
regPCrec = nwest(y,x,0); %% Full BC <- PC


retsHRec = alnrtsim_pf(2:end) .* rec_sim_02(2:end);     %% Excess Returns Recession
retsHExp = alnrtsim_pf(2:end) .* (1-rec_sim_02(2:end)); %% Exceess Returns Expansions

PDRegHRec = rec_sim_02(1:end-h,:) .* PD_regress(1:end-h,1);
PDRegHExp =  (1 - rec_sim_02(1:end-h,:)) .* PD_regress(1:end-h,1);
PCRegHRec = rec_sim_02(1:end-h,:) .* PD_regress(1:end-h,1);
PCRegHExp =  (1 - rec_sim_02(1:end-h,:)) .* PC_regress(1:end-h,1);

a = [retsHRec, PDRegHRec];
a = a(all(a,2),:);
ExcRetsRec = a(:,1);                   %% <- Excess Returns Recessions only
PDrecHR = [ones(size(a,1), 1) a(:,2)]; %% <- PD recession
regPDrec1 = nwest(ExcRetsRec,PDrecHR,0); 

a = [retsHExp, PDRegHExp];
a = a(all(a,2),:);
ExRetsExp = a(:,1);                     %% <- Excess Returns Expansions onlyu
PDexpHR   = [ones(size(a,1),1) a(:,2)]; %% <- PD Expansion
regPDexp1 = nwest(ExRetsExp,PDexpHR,0);

a = [retsHExp, PCRegHExp];
a = a(all(a,2),:);
ExRetsExp = a(:,1);
PCExpHR   = [ones(size(a,1),1) a(:,2)];
regPCexp1 = nwest(ExRetsExp,PCExpHR,0);

a = [retsHRec, PCRegHRec];
a = a(all(a,2),:);
ExRecRets = a(:,1);
PCrecHR   = [ones(size(a,1),1) a(:,2)];
regPCrec1 = nwest(ExRecRets,PCrecHR,0);

regs2 = [regPCrec regPDrec regPCrec1 regPCexp1 regPDrec1 regPDexp1];
RegressionTable1;