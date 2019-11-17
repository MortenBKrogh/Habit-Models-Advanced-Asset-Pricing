function [sg] = mkgrids(szgrid,flag)
% Will build s grid efficiently %
% ----------------------------------------------------------------------- %
global s_max S_max s_bar
if flag == 0
    sg = linspace(0,S_max,szgrid);
    aux = [(sg(szgrid)-.01) (sg(szgrid)-.02) (sg(szgrid)-.03) (sg(szgrid)-.04)];
    sg = cat(2,sg,aux);
    sg = sort(sg);
    sg = log(sg(2:size(sg,2)))'; 
    
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
    idens=log([.0005 0.0015 .0025 .0035 .0045])';
    sg=cat(1,sg,idens);
    sg = sort(sg); % Sorting the values present in the grid
end
% Grid 3 da Wachter (2005)
if flag == 1
    sg = linspace(0,S_max,szgrid);
    sg = log(sg(2:size(sg,2)))';
    if max(sg == s_bar) == 0
        sg = cat(1,sg,s_bar);
        sg = sort(sg);
    end
    if max(sg == s_max) == 0
        sg = cat(1,sg,s_max);
        sg = sort(sg);
    end
    u=min(sg);
    aux = linspace(-300,u,200)';
    sg = cat(1,sg,aux);
    sg = sort(sg);
end
end