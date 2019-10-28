function [ anss, x, w ] = GaussLegendre(f,a,b,n,tol,varargin)
% GAUSSLEGENDRE(f,a,b,n,tol) Fast and precise Gauss-Legendre quadrature.
%
% Approximate definite integral of a function f(x) on the interval [a,b]
% using n-points high precision Gauss-Legendre Quadrature.
%
% Abscissas and weights are calculated with prescribed tolerance or used
% pre-calculated with high precision.
%
% Example 1:    >>GaussLegendre(@cos,-pi/2,pi/2,1024)
%                 ans =
%                       2.000000000000000
%
% Example 2:    >>f=inline('cos(x)');
%               >>GaussLegendre(f,-pi/2,pi/2,1024)
%                 ans =
%                       2.000000000000000


% Prepare function for sending to MEX
f = fcnchk(f,'vectorized');

% Use default number of nodes and tolerance
if nargin < 4, n = 256;
end
if nargin < 5, tol = eps * 1e+3;
end

% Calculate integral by ultra-fast native compiled
% code in MEX

if nargout <= 1, anss = quadlab(20,f,a,b,n,tol);
else             [anss,x,w] = quadlab(20,f,a,b,n,tol);
end

end

