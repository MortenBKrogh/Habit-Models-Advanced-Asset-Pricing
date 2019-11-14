% This function finds a value of s_bar, which makes the percentage of the
% recession state of the economy match the empirical observed recession
% state.
%
% The function returns the difference between the empirical observed times
% of the economy being in a state of recession and the same from our
% simulation.
%
% s_bar is the steady state value of the surplus consumption ratio
%
% rec_emp is the empirical observered percentage of the economy being in a
% state of recession.
%
% astsim_pf is something we generated...
function [rec_percentage] = s_bar_match(s_bar, astsim_pf)

for i = 1:length(astsim_pf)
   
    if astsim_pf(i) < s_bar
        rec(i) = 1;
    else 
        rec(i) = 0;
        
    end
    
end

rec_percentage = sum(rec(:)==1) / length(rec) ;
end