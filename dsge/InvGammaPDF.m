function F=InvGammaPDF(x,alpha,beta);
% Luca Benati
% European Central Bank
% Monetary Policy Strategy Division
% January 2007
%
% PDF of Inverse Gamma distribution, based on Evans, Hastings, Peacock,
% Statistical Distributions, III edition
if x==0
    F=0;
else
    F=((beta^alpha)/gamma(alpha))*x^(-alpha-1)*exp(-beta/x);
end