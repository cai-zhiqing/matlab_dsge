function [alpha,beta]=GetInvGammaParameters(modeIG,stdIG)
% This function finds the parameter alpha and beta for the Inverse Gamma
% distribution, such that the mode and the standard deviation of the
% distribution are equal to modeG and stdG, respectively
% Mode: beta/(alpha-1) if alpha>1
% Variance: beta^2/((alpha-1)^2*(alpha-2))
varIG=stdIG^2;
Ratio=modeIG^2/varIG;
R=roots([1 -(4+Ratio) (5-2*Ratio) -(2+Ratio)]');
alpha=max(R);
beta=modeIG*(alpha+1);