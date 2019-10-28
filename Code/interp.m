function [fofs indx] = interp(sv,x,fx)

% Procedimento de interpolação da distribuição de s                       %
%                                                                         %
% sv -> vetor de valores quaisquer onde se quer gerar f(s).               %
% x -> vetor de pontos do grid de s. Deve ser monotônico.                 %
% fx -> vetor de valores atuais de f(x) no grid.                          %
%                                                                         %
% Acha a inclinação e o intercepto locais para ser usado sobre o grid de  %
% log(S), de tal forma que retorne f(x) = a + b*x                         %
% ----------------------------------------------------------------------- %

if isempty(min(find(fx == 0))) == 0
    
    fofs = 0*sv'; else
    T= size(x,1);
    
    % Algumas garantias sobre os inputs
    
    if x(2) < x(1)
        disp('O grid deve ser monótono e crescente');
    end
    
    if size(sv,2) > 1 && size(sv,1) == 1 sv=sv';
        chk = 1;
        
    elseif size(sv,2) > 1 && size(sv,1) > 1
        disp('ERRO: Vetor sv não pode ser uma matriz');
        
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
    
    % Estamos interessados em achar o menor índice do vetor (sv(i)-x) de forma
    % a obter um vetor indx que contenha estes índices e possamos organizar
    % monotonicamente os vetores de slope e intercepto.
    
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