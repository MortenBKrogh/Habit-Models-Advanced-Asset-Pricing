%% Regressions 2
rec_sim_02 = zeros(size(astsim,1),1);
for i = 1:size(astsim,1)
    if astsim(i) < log(0.02)
        rec_sim_02(i) = 1;
    else
        rec_sim_02(i) = 0;
    end 
end
if annual == 0
rec_sim_02 = rec_sim_02(2:end);
end
%%
y   = rets(1+h:end,1);
x   = [ones(length(rets(1:end-h,:)), 1),  ...         
    rec_sim_02(1:end-h,:) .* PC_regress(1:end-h,1), ...
    (1-rec_sim_02(1:end-h,:)) .* PC_regress(1:end-h,1)];
regPCrec2 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...            
    rec_sim_02(1:end-h,:) .* PD_regress(1:end-h,1), ...  
    (1-rec_sim_02(1:end-h,:)) .* PD_regress(1:end-h,1)]; 
regPDrec2 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...        
    (1-rec_sim_02(1:end-h,:)) .* PC_regress(1:end-h,1)];
regPCexp3 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
     rec_sim_02(1:end-h,:) .* PC_regress(1:end-h,1)];
regPCrec3 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
     rec_sim_02(1:end-h,:) .* PD_regress(1:end-h,1)];
regPDrec3 = nwest(y,x,0);

x   = [ones(length(rets(1:end-h,:)), 1),  ...         
     (1 - rec_sim_02(1:end-h,:)) .* PD_regress(1:end-h,1)];
regPDexp3 = nwest(y,x,0);

regs2 = [regPCrec2 regPDrec2 regPCrec3 regPCexp3 regPDrec3 regPDexp3];
RegressionTable1;