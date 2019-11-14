function [inside] = interd(v)
% ------------------------------------------------------------------------%
% Numerical integration expected returns P/D                              %
% ------------------------------------------------------------------------%
inside = internorm(v) .* erdinsd(v);
end