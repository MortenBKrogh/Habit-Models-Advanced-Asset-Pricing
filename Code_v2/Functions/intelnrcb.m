function out = intercb(v)
% Integrating the expected returns of public securities or                % 
% the termstructure                                           % 

global sig

out = pdf('norm',v,0,sig).*log(ercbin(v));

end