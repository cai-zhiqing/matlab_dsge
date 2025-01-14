function [BLower,BUpper]=GetLimitsSmetsWouters(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigmac_Lower = eps; 
Sigmal_Lower = eps; 
Lambda_Lower = eps; 
Phi_Lower = eps; 
ip_Lower = eps; 
iw_Lower = eps; 
Epsp_Lower = eps; 
Epsw_Lower = eps; 
Mu_p_Lower = eps; 
Mu_w_Lower = eps; 
PhiPi_Lower=1;
PhiY_Lower=0; % 0.08;
Rho_Lower=0; % 0.84;
Rho_a_Lower = eps; 
Rho_b_Lower = eps; 
Rho_g_Lower = eps; 
Rho_i_Lower = eps; 
Rho_r_Lower = eps; 
Rho_p_Lower = eps; 
Rho_w_Lower = eps; 
Sigma_a_Lower = eps; 
Sigma_b_Lower = eps; 
Sigma_g_Lower = eps; 
Sigma_i_Lower = eps; 
Sigma_r_Lower = eps; 
Sigma_p_Lower = eps; 
Sigma_w_Lower = eps; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigmac_Upper=2.5; % 1.47; %
Sigmal_Upper=3.5; % 2.30; % Elasticity of labor supply with espect to the real wage
Lambda_Upper=1; % 0.68; % Habit parameter in consumption
Phi_Upper=1; % 0.69;
ip_Upper=1-eps; % 0.21; % Indexation to past inflation
iw_Upper=1-eps; % 0.46; % Indexation to past inflation
Epsp_Upper=1-eps; % 0.73; % Price stickiness
Epsw_Upper=1-eps; % 0.74; % Wage stickiness
Mu_p_Upper=1-eps; % 0.59;
Mu_w_Upper=1-eps; % 0.62;
PhiPi_Upper=2.5; % 1.77;
PhiY_Upper=1; % 0.08;
Rho_Upper=1-eps; % 0.84;
Rho_a_Upper=1-eps; % 0.94;
Rho_b_Upper=1-eps; % 0.14;
Rho_g_Upper=1-eps; % 0.96;
Rho_i_Upper=1-eps; % 0.64;
Rho_r_Upper=1-eps; % 0.29;
Rho_p_Upper=1-eps; % 0.74;
Rho_w_Upper=1-eps; % 0.82;
Sigma_a_Upper=2; % 0.35;
Sigma_b_Upper=2; % 0.18;
Sigma_g_Upper=2; % 0.41;
Sigma_i_Upper=2; % 0.39;
Sigma_r_Upper=2; % 0.12;
Sigma_p_Upper=2; % 0.11;
Sigma_w_Upper=2; % 0.21;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLower = [Sigmac_Lower Sigmal_Lower Lambda_Lower Phi_Lower ip_Lower iw_Lower Epsp_Lower Epsw_Lower Mu_p_Lower Mu_w_Lower PhiPi_Lower PhiY_Lower Rho_Lower Rho_a_Lower Rho_b_Lower Rho_g_Lower Rho_i_Lower Rho_r_Lower Rho_p_Lower Rho_w_Lower Sigma_a_Lower Sigma_b_Lower Sigma_g_Lower Sigma_i_Lower Sigma_r_Lower Sigma_p_Lower Sigma_w_Lower]';
BUpper = [Sigmac_Upper Sigmal_Upper Lambda_Upper Phi_Upper ip_Upper iw_Upper Epsp_Upper Epsw_Upper Mu_p_Upper Mu_w_Upper PhiPi_Upper PhiY_Upper Rho_Upper Rho_a_Upper Rho_b_Upper Rho_g_Upper Rho_i_Upper Rho_r_Upper Rho_p_Upper Rho_w_Upper Sigma_a_Upper Sigma_b_Upper Sigma_g_Upper Sigma_i_Upper Sigma_r_Upper Sigma_p_Upper Sigma_w_Upper]';



    



    



