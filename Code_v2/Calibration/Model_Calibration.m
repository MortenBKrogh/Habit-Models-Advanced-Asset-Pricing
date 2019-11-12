function [Coefficients] = Model_Calibration
Coefficients = struct;
%% Risk free rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The Risk free rate is calculated from annual returns on  %%%
%%% 90 day T-bills less the inflation rate. Using CRSP index %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RR = importdata('Bonds_INF.csv');
Bond = RR.data(:,2);
Infl = RR.data(:,4);
Rbond = Bond - Infl;
RFR = log(Rbond + 1); % Risk free rate calib; 
meanRFR = mean(RFR);
Coefficients.rf = meanRFR;
clearvars -except Coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% the mean of the risk free rate 1950-2019 is found to be .0109 %%%
%%% that is slightly above 1%                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Persistence coefficient
RetAM = importdata('All_returns_market_Monthly.csv');
R =  1+RetAM.data(:,2); 
Rx =  1+RetAM.data(:,3); 
pd = log(R./Rx-1);
pd_t  = pd(1:end-1,:); % t
pd_t1 = pd(2:end,:); % t+1
AR1 = pd_t1\[ones(size(pd_t,1),1) pd_t];
Coefficients.Phi = AR1(2)^12;
%clearvars -except Coefficients R Rx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Autocorrelation of price/dividend ratio is found to be .9008  %%%
%%% that is slightly above .87                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Consumption growth moments
%url = 'https://fred.stlouisfed.org/';
%c = fred(url);
%startdate = '01/01/1950';
%enddate = floor(now);
%cons = fetch(c, 'A796RX0Q048SBEA',startdate,enddate); % Real Consumption non-durable goods
%consdat = cons.Data(:,2);
%writematrix(consdat);
dat = readtable('consdat.txt');
dat = table2array(dat);
dat = log(dat);
diffDat = dat(2:end) - dat(1:end-1);
g = mean(diffDat)*4;
Coefficients.g = g;
Coefficients.sigma = std(diffDat) * sqrt(4);
clearvars -except Coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Multiplying with 4 as we use quarterly data here, g is found to be %%%
%%% .0134 or 1.34% while sigma_v is found to be 1.52% (.0152)          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standard deviation of dividend growth
RetAM = importdata('All_returns_market_Quart.csv');
RetAM = RetAM.data(:,[2 3]);
Divs = (RetAM(:,1) - RetAM(:,2));
LogDivs = log(Divs);
DiffDiv = diff(LogDivs)./3;
sigma_w = std(DiffDiv)*sqrt(12);
Coefficients.sigma_w = sigma_w;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Quarterly Data on returns yields a sigma_w f .1777 or 17.77%       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
