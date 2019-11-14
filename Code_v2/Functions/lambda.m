function [y]=lambda(s)
% function Lambda
global S_bar s_bar s_max verd
if verd == 0
    if (double(s) <= s_max)
        y = (1 / S_bar).*sqrt(max(0 , 1-2.*(s-s_bar)))-1;
    else
        y = 0;
    end
elseif verd == 1
    y=(1-S_bar)/S_bar;
end
end