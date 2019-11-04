function [fofs indx] = interp(sv,x,fx)
% ------------------------------------------------------------------------%
% Interpolation procedure of s distribution                               %
% sv -> vector of any values where f (s) is to be generated.              %
% x -> vector of grid points of s. It must be monotonic.                  %
% fx -> vector of current values of f (x) in the grid.                    %
% Find local slope and intercept to use on% grid                          %
% log (S) such that returns f (x) = a + b * x                             %
% ------------------------------------------------------------------------%

if isempty(min(find(fx == 0))) == 0    
    fofs = 0*sv'; else
    T= size(x,1);
    
    if x(2) < x(1)
        disp('O grid deve ser mon?tono e crescente');
    end
    
    if size(sv,2) > 1 && size(sv,1) == 1 sv=sv';
        chk = 1;
        
    elseif size(sv,2) > 1 && size(sv,1) > 1
        disp('ERRO: Vetor sv n?o pode ser uma matriz');
        
    elseif size(sv,2) == 1 && size(sv,1) > 1
        chk = 0;
        
    elseif size(sv,2) == 1 && size(sv,1) == 1
        chk = 2;                                  % Para o caso de sv = s %
        
    end
    
    gradf = (fx(2:T)-fx(1:T-1))./(x(2:T)-x(1:T-1));
    const = cat(1,(fx(1:T-1)-gradf.*x(1:T-1)),(fx(T)-gradf(T-1)*x(T)));
    const = cat(1,fx(1)-gradf(1)*x(1),const);
    slope = cat(1,gradf(1),gradf);
    slope = cat(1,slope,gradf(T-1));
    
  % We are interested in finding the smallest vector index (sv (i) -x) so
  % we get an index vector that contains these indexes and we can arrange
  % monotonically the slope and intercept vectors.

    indx = zeros(size(sv,1),1);
    
    for i = 1 : size(sv,1)
        
        if sv(i)> x(1)
            indx(i)= max(find((sv(i)-x) > 0))+1;
        else
            indx(i)=1;
        end
    end
    
    fofs = const(indx)+ slope(indx).*sv;
    if chk == 0
        fofs = fofs';
    end
end
end