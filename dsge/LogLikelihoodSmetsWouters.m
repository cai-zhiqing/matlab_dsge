function [MinusLogLikelihood,L]=LogLikelihoodSmetsWouters(x,Y,T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Beta=0.99; % Discount factor applied by households
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
PhiPi=x(11); % r_pi in Smets paper equation 14
PhiY=x(12);  % r_y in Smets paper equation 14
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gammabar=0.44;
ep=10; % Curvature of Kimball goods market aggregator
ew=10; % Curvature of Kimball labor market aggregator
Phip=1.6; % 1 + share of fixed costs in production
Phiw=1.5; % 1 + steady-state labor market markup
Alpha=0.21; % Share of capital in production
Kappay=1.7; % Steady-state capital-output ratio, based on long-run data, from D'Adda and Scorcu
PhiDY=0.16;
Rho_ga=0.39;
Delta=0.025; % Depreciation rate of capital
Csi=6.23; % Steady-state elasticity of capital adjustment cost function
Gamma=1+Gammabar/100; % Steady-state growth rate
Rkstar=Gamma^Sigmac/Beta-(1-Delta); % Steady-state rental rate if capital
zy=Rkstar*Kappay;
WhstarTimeLstarOverCstar=1.1111;
gy=0.18; % Steady-state spending-output ratio
iy=(Gamma-1+Delta)*Kappay; % Steady-state investment-output ratio
cy=1-gy-iy; % Steady-state share of consumption in output
%
c1=(Lambda/Gamma)/(1+(Lambda/Gamma));
c2=((Sigmac-1)*(WhstarTimeLstarOverCstar))/(Sigmac*(1+(Lambda/Gamma)));
c3=(1-Lambda/Gamma)/(Sigmac*(1+(Lambda/Gamma)));
%
i1=1/(1+Beta*Gamma^(1-Sigmac));
i2=1/((1+Beta*Gamma^(1-Sigmac))*Gamma^2*Csi);
%
q1=(1-Delta)*Beta*Gamma^(-Sigmac);
%
z1=(1-Phi)/Phi;
%
k1=(1-Delta)/Gamma;
k2=(1-(1-Delta)/Gamma)*(1+Beta*Gamma^(1-Sigmac))*Gamma^2*Csi;
%
pi1=ip/(1+Beta*Gamma^(1-Sigmac)*ip);
pi2=Beta*Gamma^(1-Sigmac)/(1+Beta*Gamma^(1-Sigmac)*ip);
pi3=(1/(1+Beta*Gamma^(1-Sigmac)*ip))*(((1-Beta*Gamma^(1-Sigmac)*Epsp)*(1-Epsp))/(Epsp*(1+ep*(Phip-1))));
%
w1=1/(1+Beta*Gamma^(1-Sigmac));
w2=(1+Beta*Gamma^(1-Sigmac)*iw)/(1+Beta*Gamma^(1-Sigmac));
w3=iw/(1+Beta*Gamma^(1-Sigmac));
w4=(1/(1+Beta*Gamma^(1-Sigmac)))*(((1-Beta*Gamma^(1-Sigmac)*Epsw)*(1-Epsw))/(Epsw*(1+ew*(Phiw-1))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                The matrices in the Sims canonical form:
%
N=29;
%
G0=zeros(N);
G0(14:20,1:7)=eye(7); G1(21:N,21:N)=eye(9);
G0(1,:)=[-cy zeros(1,2) -iy 0 -zy*z1 0 1 zeros(1,12) -1 zeros(1,8)];
G0(2,:)=[1 -c2 zeros(1,10) c3 c1-1 c2 -c3 zeros(1,5) c3 zeros(1,7)];
G0(3,:)=[zeros(1,3) 1 -i2 zeros(1,11) i1-1 zeros(1,5) -1 zeros(1,6)];
G0(4,:)=[zeros(1,4) 1 zeros(1,7) 1 zeros(1,2) -1 0 -q1 (q1-1) zeros(1,2) 1 zeros(1,7)];
G0(5,:)=[0 (Alpha-1)*Phip zeros(1,5) 1 -Alpha*Phip zeros(1,14) -Phip zeros(1,5)];
G0(6,:)=[zeros(1,5) -z1 zeros(1,2) 1 zeros(1,20)];
G0(7,:)=[zeros(1,3) k1-1 zeros(1,5) 1 zeros(1,12) -k2 zeros(1,6)];
G0(8,:)=[0 Alpha zeros(1,4) 1 0 -Alpha 0 1 zeros(1,12) -1 zeros(1,5)];
G0(9,:)=[zeros(1,2) 1 zeros(1,7) pi3 zeros(1,4) -pi2 zeros(1,8) -1 zeros(1,4)];
G0(10,:)=[0 -1 zeros(1,3) 1 -1 0 0 1 zeros(1,19)];
G0(11,:)=[1/(1-Lambda/Gamma) Sigmal zeros(1,4) -1 zeros(1,4) 1 zeros(1,17)];
G0(12,:)=[zeros(1,2) w2 zeros(1,3) 1 zeros(1,4) w4 zeros(1,3) (w1-1) zeros(1,3) (w1-1) zeros(1,5) -1 zeros(1,3)];
G0(13,:)=[zeros(1,2) (Rho-1)*PhiPi zeros(1,4) -((1-Rho)*PhiY+PhiDY) zeros(1,4) 1 zeros(1,13) -1 zeros(1,2)];
G0(21:N,21:N)=eye(9);
%
G1=zeros(N); G1(2,1)=c1; G1(3,4)=i1; G1(6,10)=1; G1(7,10)=k1; G1(9,3)=pi1; G1(12,3)=w3; G1(12,7)=w1; G1(11,1)=(Lambda/Gamma)/(1-(Lambda/Gamma));
G1(13,8)=-PhiDY; G1(13,13)=Rho; G1(14:20,14:20)=eye(7); G1(21:27,21:27)=diag([Rho_g Rho_b Rho_i Rho_a Rho_p Rho_w Rho_r]'); G1(25:26,28:29)=diag([-Mu_p -Mu_w]');
%
CSI=[zeros(20,7); eye(7); [zeros(2,4) eye(2) zeros(2,1)]];
CSI(21,4)=Rho_ga;
%
PI=[zeros(13,7); eye(7); zeros(9,7)];
%
% The QZ decomposition of the matrix pencil (G0-mu*G1)
%
[L,O,Q,Z]=qz(G0,G1);
[L,O,Q,Z]=qzdiv(1.0001,L,O,Q,Z); % Reorder the decomposition so that all the explosive eigenvalues are in the lower right corner.
lambda=diag(O)./diag(L); % The eigenvalues
INST=sum(abs(lambda)>=1); % The number of unstable eigenvalues
%
CSIstar=Q*CSI;
PIstar=Q*PI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [8] y(t)
% [1] c(t)
% [4] i(t)
% [7] w(t)
% [2] l(t)
% [3] inf(t)
% [13] r(t)
NumberOfObservables=7;
N=size(G0,1);
%
H=zeros(NumberOfObservables,N);
H(1,8)=1;
H(2,1)=1;
H(3,4)=1;
H(4,7)=1;
H(5,2)=1;
H(6,3)=1;
H(7,13)=1;
H=H';
SelectorMatrix=H;
%
VAR=diag([Sigma_g Sigma_b Sigma_i Sigma_a Sigma_p Sigma_w Sigma_r]'.^2);
%
PIstar=Q*PI;
PIstarx=PIstar(N+1-INST:N,:);
PIstarc=PIstar(1:N-INST,:);
CSIstar=Q*CSI;
CSIstarc=CSIstar(1:N-INST,:);
CSIstarx=CSIstar(N+1-INST:N,:);
L11=L(1:N-INST,1:N-INST);
O11=O(1:N-INST,1:N-INST);
Z1=Z(:,1:N-INST);
%
if INST==7 % Here we have determinacy
    if max(max(abs(imag(Z1*(L11\(CSIstarc-PIstarc*(PIstarx\CSIstarx)))))))>1.0e-10
        MinusLogLikelihood=1.0e+100;
        L=-9999;
        return
    end
    ImpactMatrix=real(Z1*(L11\(CSIstarc-PIstarc*(PIstarx\CSIstarx))));
else
    MinusLogLikelihood=1.0e+100;
    L=-9999;
    return
end
Q=real(ImpactMatrix*VAR*ImpactMatrix');
F=real(Z1*(L11\O11)*Z1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This is the Kalman filtering recursion to compute the log-likelihood:
%
% Initial values for state vector:
s0=zeros(size(F,1),1);
% Initial values for precision matrix of state vector:
P0=eye(size(Q))*5^2;
% This is s(t|t-1):
s=zeros(size(s0,1),T+1);
% This is s(t|t):
st=s;
% Initialize the first observation of the state vector:
s(:,1)=s0;
P=P0;
if det(H'*P*H)==0
    MinusLogLikelihood=1.0e+100;
    L=-9999;
    return
end
st(:,1)=s0;
Pt=P0;
% Pre-allocate the space for the the log-likelihood:
L=zeros(T,1);
% The Kalman filter recursion:     
tt=1;
while tt<=T
    % The one-step ahead forecast error for the observed variables:
    uhat=Y(tt,:)'-H'*s(:,tt);
    % This is s(t|t):
    st(:,tt+1)=s(:,tt)+P*H*MyInverse(H'*P*H)*uhat;
    % This is its precision matrix, P(t|t):
    Pt=P-P*H*(MyInverse(H'*P*H))*H'*P;
    % The log-likelihood for observation t:
    L(tt)=-0.5*log(det(H'*P*H))-0.5*(uhat'*(MyInverse(H'*P*H))*uhat);
    % This is s(t+1|t):
    s(:,tt+1)=F*s(:,tt)+F*P*H*(MyInverse(H'*P*H))*uhat;
    % This is the precision matrix of s(t+1|t), P(t+1|t):
    P=Q+F*P*F'-F*P*H*(MyInverse(H'*P*H))*H'*P*F';
    if det(H'*P*H)==0
        MinusLogLikelihood=1.0e+100;
        return
    end
    tt=tt+1;
end
if sum(abs(imag(L))>0)>0
    MinusLogLikelihood=1.0e+100;
else
    MinusLogLikelihood=-sum(L);
end
%





