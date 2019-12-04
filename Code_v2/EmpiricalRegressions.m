addpath('Data');
addpath('Calibration');
addpath('Functions');
clear
%%
datfile = importdata('All_returns_market_Monthly.csv');
%datfile = importdata('All_returns_market_Annual.csv');
cal = datfile.data(:,1); % Dates
r  = datfile.data(:,2); % Returns
rx  = datfile.data(:,3); % Returns less dividinds
%%
dp =  (1+r)./(1+rx) - 1; % Div. Yield
dd = [NaN ; dp(2:end)./dp(1:end-1).*(1+rx(2:end))-1];  % Div. Growth
%%
RR = importdata('Rets30dayTBILL.csv');
Bond = RR.data(:,2); % Returns on 90 day T-bill
Infl = RR.data(:,3); % Inflation rate measured as CPI
RFR = Bond - Infl; % Real Risk free rate
re = r-RFR; % Excess Returns
%%
RecData = importdata('USREC_Period.csv');
Rec = RecData.data(:,1);
%%
T = length(dp);
dpRec = dp .* Rec;
reRec = re .* Rec;

a = [dpRec(1:end-1), reRec(2:end)];
a = a(all(a,2),:);

rhv = [ones(length(a),1), a(:,1)]; 
lhv = a(:,2);
Rec_Reg = nwest(lhv,rhv,1)
%%
T = length(dp);
dpExp = dp .* (1-Rec);
reExp = re .* (1-Rec);

a = [dpExp(1:end-1), reExp(2:end)];
a = a(all(a,2),:);

rhv = [ones(length(a),1), a(:,1)]; 
lhv = a(:,2);
Rec_Reg = nwest(lhv,rhv,1)