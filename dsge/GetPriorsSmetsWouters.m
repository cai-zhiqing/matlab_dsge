function [A,B]=GetPriorsSmetsWouters(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mode_Sigmac = 1.39;
Std_Sigmac=0.1;
[A_Sigmac,B_Sigmac]=GetGammaParameters(Mode_Sigmac,Std_Sigmac);
%
Mode_Sigmal = 1.92;
Std_Sigmal=0.1;
[A_Sigmal,B_Sigmal]=GetGammaParameters(Mode_Sigmal,Std_Sigmal);
%
Mode_Lambda = 0.71;
Std_Lambda = 0.1;
[A_Lambda,B_Lambda]=GetBetaParameters(Mode_Lambda,Std_Lambda);
%
Mode_Phi = 0.54;
Std_Phi = 0.1;
[A_Phi,B_Phi]=GetBetaParameters(Mode_Phi,Std_Phi);
%
Mode_ip = 0.22;
Std_ip = 0.1;
[A_ip,B_ip]=GetBetaParameters(Mode_ip,Std_ip);
%
Mode_iw = 0.59;
Std_iw = 0.1;
[A_iw,B_iw]=GetBetaParameters(Mode_iw,Std_iw);
%
Mode_Epsp = 0.65;
Std_Epsp = 0.1;
[A_Epsp,B_Epsp]=GetBetaParameters(Mode_Epsp,Std_Epsp);
%
Mode_Epsw = 0.73;
Std_Epsw = 0.1;
[A_Epsw,B_Epsw]=GetBetaParameters(Mode_Epsw,Std_Epsw);
%
Mode_Mu_p = 0.74;
Std_Mu_p = 0.1;
[A_Mu_p,B_Mu_p]=GetBetaParameters(Mode_Mu_p,Std_Mu_p);
%
Mode_Mu_w = 0.88;
Std_Mu_w = 0.1;
[A_Mu_w,B_Mu_w]=GetBetaParameters(Mode_Mu_w,Std_Mu_w);
%
Mode_PhiPi = 1.5;
Std_PhiPi=0.1;
[A_PhiPi,B_PhiPi]=GetGammaParameters(Mode_PhiPi,Std_PhiPi);
%
Mode_PhiY = 0.5;
Std_PhiY=0.1;
[A_PhiY,B_PhiY]=GetGammaParameters(Mode_PhiY,Std_PhiY);
%
Mode_Rho = 0.95;
Std_Rho = 0.1;
[A_Rho,B_Rho]=GetBetaParameters(Mode_Rho,Std_Rho);
%
Mode_Rho_a = 0.95;
Std_Rho_a = 0.02;
[A_Rho_a,B_Rho_a]=GetBetaParameters(Mode_Rho_a,Std_Rho_a);
%
Mode_Rho_b = 0.18;
Std_Rho_b = 0.1;
[A_Rho_b,B_Rho_b]=GetBetaParameters(Mode_Rho_b,Std_Rho_b);
%
Mode_Rho_g = 0.97;
Std_Rho_g = 0.01;
[A_Rho_g,B_Rho_g]=GetBetaParameters(Mode_Rho_g,Std_Rho_g);
%
Mode_Rho_i = 0.71;
Std_Rho_i = 0.1;
[A_Rho_i,B_Rho_i]=GetBetaParameters(Mode_Rho_i,Std_Rho_i);
%
Mode_Rho_r = 0.12;
Std_Rho_r = 0.05;
[A_Rho_r,B_Rho_r]=GetBetaParameters(Mode_Rho_r,Std_Rho_r);
%
Mode_Rho_p = 0.9;
Std_Rho_p = 0.05;
[A_Rho_p,B_Rho_p]=GetBetaParameters(Mode_Rho_p,Std_Rho_p);
%
Mode_Rho_w = 0.97;
Std_Rho_w = 0.01;
[A_Rho_w,B_Rho_w]=GetBetaParameters(Mode_Rho_w,Std_Rho_w);
%
Mode_Sigma_a = 0.45;
Std_Sigma_a = 0.1;
[A_Sigma_a,B_Sigma_a]=GetInvGammaParameters(Mode_Sigma_a,Std_Sigma_a);
%
Mode_Sigma_b = 0.24;
Std_Sigma_b = 0.1;
[A_Sigma_b,B_Sigma_b]=GetInvGammaParameters(Mode_Sigma_b,Std_Sigma_b);
%
Mode_Sigma_g = 0.52;
Std_Sigma_g = 0.1;
[A_Sigma_g,B_Sigma_g]=GetInvGammaParameters(Mode_Sigma_g,Std_Sigma_g);
%
Mode_Sigma_i = 0.45;
Std_Sigma_i = 0.1;
[A_Sigma_i,B_Sigma_i]=GetInvGammaParameters(Mode_Sigma_i,Std_Sigma_i);
%
Mode_Sigma_r = 0.24;
Std_Sigma_r = 0.1;
[A_Sigma_r,B_Sigma_r]=GetInvGammaParameters(Mode_Sigma_r,Std_Sigma_r);
%
Mode_Sigma_p = 0.12;
Std_Sigma_p = 0.1;
[A_Sigma_p,B_Sigma_p]=GetInvGammaParameters(Mode_Sigma_p,Std_Sigma_p);
%
Mode_Sigma_w = 0.24;
Std_Sigma_w = 0.1;
[A_Sigma_w,B_Sigma_w]=GetInvGammaParameters(Mode_Sigma_w,Std_Sigma_w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Encoding all the priors in two vectors A and B:
%
A=[A_Sigmac A_Sigmal A_Lambda A_Phi A_ip A_iw A_Epsp A_Epsw A_Mu_p A_Mu_w A_PhiPi A_PhiY A_Rho A_Rho_a A_Rho_b A_Rho_g A_Rho_i A_Rho_r A_Rho_p A_Rho_w A_Sigma_a A_Sigma_b A_Sigma_g A_Sigma_i A_Sigma_r A_Sigma_p A_Sigma_w]';
B=[B_Sigmac B_Sigmal B_Lambda B_Phi B_ip B_iw B_Epsp B_Epsw B_Mu_p B_Mu_w B_PhiPi B_PhiY B_Rho B_Rho_a B_Rho_b B_Rho_g B_Rho_i B_Rho_r B_Rho_p B_Rho_w B_Sigma_a B_Sigma_b B_Sigma_g B_Sigma_i B_Sigma_r B_Sigma_p B_Sigma_w]';







