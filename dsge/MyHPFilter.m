function [yT,yC]=MyHPFilter(y,lambda)
%
% Luca Benati
% Monetary Policy Strategy Division
% European Central Bank
%
% function [yT,yC]=MyHPFilter(y,lambda)
% This program HP-filters a previously defined time series y.
% The algorithm used in the program is based on Danthine and
% Girardin, 'European Economic Review', 1989.
%
l=length(y);
K=zeros(l-2,l);
for j=1:(l-2)
    K(j,j)=1;
    K(j,j+1)=-2;
    K(j,j+2)=1;
end
A=eye(l)+lambda*(K'*K);
yT=A\y;
yC=y-yT;
if nargout==1
    yT=yC;
end