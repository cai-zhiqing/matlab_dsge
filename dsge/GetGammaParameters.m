function [kappa,theta]=GetGammaParameters(modeG,stdG)
% This function finds the parameter kappa and theta for the Gamma
% distribution, such that the mode and the standard deviation of the
% distribution are equal to modeG and stdG, respectively
% Mode: (kappa-1)*theta
% Variance: kappa*theta^2
VarG=stdG^2;
R=roots([VarG -(2*VarG+modeG^2) VarG]');
kappa=max(R);
theta=modeG/(kappa-1);