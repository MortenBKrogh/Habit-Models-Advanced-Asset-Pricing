function [fx] = intgl1(f,a,b)
% Rotina que desenvolve o método de integração de Gauss-Legendre. %
% Nesta rotina, o vetor 'e' compõe-se de valores entre 0 e 1 (e ~ (0,1)) % e o vetor 'w' das raízes do polinômio de Legendre cujos coeficientes vem % de 'e'. Estes valores estão contidos no arquivo ewquad.txt e foram
% processados pelo programa 'gausslege.m' que não faz parte deste trabalho % por fugir muito do seu escopo.
[e w] = textread('ewquad.txt','%f %f');
e = sort(cat(1,e,-e));
w = cat(1,sort(w),w);

diff = b - a;
xc = 0.5*( (b+a)+(diff * e'));
fx = zeros(size(xc)) + f(xc);
fx = ((diff/2)* (fx*w));
end