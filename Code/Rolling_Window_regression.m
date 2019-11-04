clear
clc
%%
Dat = readtable('ts1.txt');
tsc = 12;
g=0.0189/tsc;
sig=0.015/sqrt(tsc);
rf0=0.0094/tsc;
phi=0.87^(1/tsc);
gamma=2;
B=0;
verd=0;
ann=0;
S_bar=sig*sqrt(gamma/(1-phi-B/gamma));
s_bar = log(S_bar);
s_max = s_bar + (1-S_bar^2)/2;
S_max = exp(s_max);
%% Variables
s_t = Dat.S_t;
pd_t = Dat.PCratio;
ret_t = Dat.ExPostReturns;
%% State Variables
T = size(s_t, 1);
S_state_Rec = zeros(T, 1);
S_state_Exp = zeros(T, 1);
Dummy_Rec   = zeros(T, 1);
for i=1:size(s_t, 1)
    if s_t(i) > s_bar
        S_state_Exp(i) = s_t(i);
        S_state_Rec(i) = 0;
        Dummy_Rec(i)   = 0;
    else
        S_state_Exp(i) = 0;
        S_state_Rec(i) = s_t(i);
        Dummy_Rec(i)   = 1;
    end
end
%%
Rolls = 3000;
nahead = 1;
beta  = NaN(Rolls, 3);
stats = NaN(Rolls, 4);
std = NaN(3,2,Rolls);
for i = 1:Rolls
    Mod = ret_t(nahead+i:i+4999+nahead,:);
    pred = [ones(size(Mod,1), 1), ret_t(i:i+4999), pd_t(i:i+4999) .* Dummy_Rec(i:i+4999)];
    [beta1, std1, ~, ~, stats1] = regress(Mod, pred);
    beta(i,:)  = beta1;
    std(:,:,i)   = std1;
    stats(i,:) = stats1;
end
%%
 plot(stats(:,1))
