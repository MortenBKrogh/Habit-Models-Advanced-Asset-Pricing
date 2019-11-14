function [out] = inter2d(v)

% Expected variance returns consumption claim
% ----------------------------------------------------------------------- %

global sig g rhow sig_w lnpca sg
out = internorm(v) .* erd2ind(v);
end