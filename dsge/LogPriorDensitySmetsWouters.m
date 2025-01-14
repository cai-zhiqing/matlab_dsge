function LPD=LogPriorDensitySmetsWouters(x,A,B,NoLog)
if nargin<4
    NoLog='N';
end
x=x(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigmac=x(1);
Sigmal=x(2);
Lambda=x(3);
Phi=x(4);
ip=x(5);
iw=x(6);
Epsp=x(7);
Epsw=x(8);
Mu_p=x(9);
Mu_w=x(10);
PhiPi=x(11);
PhiY=x(12);
Rho=x(13);
Rho_a=x(14);
Rho_b=x(15);
Rho_g=x(16);
Rho_i=x(17);
Rho_r=x(18);
Rho_p=x(19);
Rho_w=x(20);
Sigma_a=x(21);
Sigma_b=x(22);
Sigma_g=x(23);
Sigma_i=x(24);
Sigma_r=x(25);
Sigma_p=x(26);
Sigma_w=x(27);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_Sigmac=A(1);
A_Sigmal=A(2);
A_Lambda=A(3);
A_Phi=A(4);
A_ip=A(5);
A_iw=A(6);
A_Epsp=A(7);
A_Epsw=A(8);
A_Mu_p=A(9);
A_Mu_w=A(10);
A_PhiPi=A(11);
A_PhiY=A(12);
A_Rho=A(13);
A_Rho_a=A(14);
A_Rho_b=A(15);
A_Rho_g=A(16);
A_Rho_i=A(17);
A_Rho_r=A(18);
A_Rho_p=A(19);
A_Rho_w=A(20);
A_Sigma_a=A(21);
A_Sigma_b=A(22);
A_Sigma_g=A(23);
A_Sigma_i=A(24);
A_Sigma_r=A(25);
A_Sigma_p=A(26);
A_Sigma_w=A(27);
%
B_Sigmac=B(1);
B_Sigmal=B(2);
B_Lambda=B(3);
B_Phi=B(4);
B_ip=B(5);
B_iw=B(6);
B_Epsp=B(7);
B_Epsw=B(8);
B_Mu_p=B(9);
B_Mu_w=B(10);
B_PhiPi=B(11);
B_PhiY=B(12);
B_Rho=B(13);
B_Rho_a=B(14);
B_Rho_b=B(15);
B_Rho_g=B(16);
B_Rho_i=B(17);
B_Rho_r=B(18);
B_Rho_p=B(19);
B_Rho_w=B(20);
B_Sigma_a=B(21);
B_Sigma_b=B(22);
B_Sigma_g=B(23);
B_Sigma_i=B(24);
B_Sigma_r=B(25);
B_Sigma_p=B(26);
B_Sigma_w=B(27);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_Sigmac=Gammapdf(Sigmac,A_Sigmac,B_Sigmac);
p_Sigmal=Gammapdf(Sigmal,A_Sigmal,B_Sigmal);
p_Lambda=betapdf(Lambda,A_Lambda,B_Lambda);
p_Phi=betapdf(Phi,A_Phi,B_Phi);
p_ip=betapdf(ip,A_ip,B_ip);
p_iw=betapdf(iw,A_iw,B_iw);
p_Epsp=betapdf(Epsp,A_Epsp,B_Epsp);
p_Epsw=betapdf(Epsw,A_Epsw,B_Epsw);
p_Mu_p=betapdf(Mu_p,A_Mu_p,B_Mu_p);
p_Mu_w=betapdf(Mu_w,A_Mu_w,B_Mu_w);
p_PhiPi=Gammapdf(PhiPi,A_PhiPi,B_PhiPi);
p_PhiY=Gammapdf(PhiY,A_PhiY,B_PhiY);
p_Rho=betapdf(Rho,A_Rho,B_Rho);
p_Rho_a=betapdf(Rho_a,A_Rho_a,B_Rho_a);
p_Rho_b=betapdf(Rho_b,A_Rho_b,B_Rho_b);
p_Rho_g=betapdf(Rho_g,A_Rho_g,B_Rho_g);
p_Rho_i=betapdf(Rho_i,A_Rho_i,B_Rho_i);
p_Rho_r=betapdf(Rho_r,A_Rho_r,B_Rho_r);
p_Rho_p=betapdf(Rho_p,A_Rho_p,B_Rho_p);
p_Rho_w=betapdf(Rho_w,A_Rho_w,B_Rho_w);
p_Sigma_a=InvGammaPDF(Sigma_a,A_Sigma_a,B_Sigma_a);
p_Sigma_b=InvGammaPDF(Sigma_b,A_Sigma_b,B_Sigma_b);
p_Sigma_g=InvGammaPDF(Sigma_g,A_Sigma_g,B_Sigma_g);
p_Sigma_i=InvGammaPDF(Sigma_i,A_Sigma_i,B_Sigma_i);
p_Sigma_r=InvGammaPDF(Sigma_r,A_Sigma_r,B_Sigma_r);
p_Sigma_p=InvGammaPDF(Sigma_p,A_Sigma_p,B_Sigma_p);
p_Sigma_w=InvGammaPDF(Sigma_w,A_Sigma_w,B_Sigma_w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LPD=p_Sigmac*p_Sigmal*p_Lambda*p_Phi*p_ip*p_iw*p_Epsp*p_Epsw*p_Mu_p*p_Mu_w*p_PhiPi*p_PhiY*p_Rho*p_Rho_a*p_Rho_b*p_Rho_g*p_Rho_i*p_Rho_r*p_Rho_p*p_Rho_w*p_Sigma_a*p_Sigma_b*p_Sigma_g*p_Sigma_i*p_Sigma_r*p_Sigma_p*p_Sigma_w;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if NoLog=='N' 
    % This is the log of the prior density:
    LPD=-log(LPD);
else
    % This is just the prior density itself:
    LPD=LPD;
end


    



        
        


    


