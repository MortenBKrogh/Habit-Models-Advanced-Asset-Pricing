function xs=selif(x,t)
% SELIF 
% SELIF(x,t) selects the elements of x for which t=1
% x: NxK matrix, t: Nx1 matrix of 0's and 1's
nt=repmat(t,1,size(x',1));
nx=x.*nt;
xs=nx(nx~=0);
col=size(x',1);
xs=reshape(xs,size(xs,1)/col,col);
