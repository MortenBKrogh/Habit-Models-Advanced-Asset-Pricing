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
RR = importdata('Bonds_INF.csv');
Bond = RR.data(:,2); % Returns on 90 day T-bill
Infl = RR.data(:,4); % Inflation rate measured as CPI
RFR = Bond - Infl; % Real Risk free rate
re = r-RFR; % Excess Returns
%%
RecData = importdata('USREC.csv');
Rec = Rec.data(:,1);
%%
T = length(dp);
rhv = [ones(T-1,1) dp(1:T-1)]; 
lhv = re(2:end); 
nwest(lhv,rhv,0)