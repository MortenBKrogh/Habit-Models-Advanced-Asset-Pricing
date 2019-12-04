addpath('Data');
addpath('Calibration');

clear
%%
datfile = importdata('All_returns_market_Monthly.csv');
retdat  = datfile.data(:,2:3);
