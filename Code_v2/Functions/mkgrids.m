function [sg] = mkgrids(szgrid)
global s_max s_bar
a = exp(s_max)*1;
b = (exp(s_max)*1)/(szgrid+1);
    sg = 0:b:a;
    sg = log(sg(2:end))';      
    
     if max(sg == s_bar) == 0       % Making sure s_bar will be on the grid.
        sg = cat(1,sg,s_bar);
        sg = sort(sg);
    end
    if max(sg == s_max) == 0
        sg = cat(1,sg,s_max);
        sg = sort(sg);
    end
                                            
    % Put more density at the beginning of the grid to improve iteration
    % during the fixed point procedure.%
    % ----------------------------------------------------------------------- %
    idens= s_max-[0.01:0.01:0.04]';
    sg=cat(1,sg,idens);
    sg = sort(sg); % Sorting the values present in the grid
end