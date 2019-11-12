function [y]=lambda_Helper(s)
% function Lambda
global S_bar s_bar s_max verd
        y = (1 / S_bar).*sqrt(max(0 , 1-2.*(s-s_bar)))-1;
end