%% Variables:

%%% Consumption %%%
% ln(C_t) or c_t
ct = cumsum(alndctsim_pf);
% delta c_t
dct = alndctsim_pf;

%%% Surplus Consumption ratio %%%
% ln(S_t) or s_t
st = astsim_pf;
% S_t
St = exp(astsim_pf);

%%% Price/Consumption ratio %%%
% ln(P_t/C_t) or p_t - c_t
pct = alnpctsim_pf;

%%% Returns %%%
% ln(R_t) or r_t
rt = alnrtsim_pf;

%  std(r_t) or Volatility
std_rt = asdlnrtsim_pf;

%%% Risk-free rate %%%
% log(R^f _t) or r^f_t
rft = alnrfsim_pf;

%%% LOM: Prices %%%
% (p_t - c_t) - (p_(t-1) - c_(t-1)) + delta c_t
pt = alnchpsim_pf;
