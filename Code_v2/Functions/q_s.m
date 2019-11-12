function [q] = q_s(s)
%Density of s continous time
global s_max
q = z_s(s)./integral(@z_s,-Inf,s_max);
end

