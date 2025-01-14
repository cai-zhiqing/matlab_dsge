%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
options = optimset;
options = optimset(options,'Display','off');
warning('off','all')
%path(path,'U:\SeminarEmpiricalMacro')
path(path,'D:\code_github\matlab\dsge')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Country='EA';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Country=='US'
    X=xlsread('D:\code_github\matlab\dsge\Dataset.xlsx','US','A41:H258');
    [T,N]=size(X);
    Time = X(2:T,1);
    ShortRate = X(:,2);
    Inflation=diff(log(X(:,3)))*400;
    OutputGap = MyHPFilter(log(X(:,4)),1600)*100;
    HPFilteredConsumption = MyHPFilter(X(:,5),1600)*100;
    HPFilteredInvestment = MyHPFilter(X(:,6),1600)*100;
    HPFilteredHours=MyHPFilter(log(X(:,7)),1600)*100;
    HPFilteredRealWage=MyHPFilter(log(X(:,8)),1600)*100;
elseif Country=='EA'
    X=xlsread('D:\code_github\matlab\dsge\Dataset.xlsx','EuroArea','A12:H211');
    [T,N]=size(X);
    Time = X(2:T,1);
    ShortRate = X(:,2);
    Inflation=diff(log(X(:,3)))*400;
    OutputGap = MyHPFilter(log(X(:,4)),1600)*100;
    HPFilteredConsumption = MyHPFilter(log(X(:,5)),1600)*100;
    HPFilteredInvestment = MyHPFilter(log(X(:,6)),1600)*100;
    HPFilteredHours=MyHPFilter(log(X(:,7)),1600)*100;
    HPFilteredRealWage=MyHPFilter(log(X(:,8)),1600)*100;
else
end
Y=demean([OutputGap(2:end) HPFilteredConsumption(2:end) HPFilteredInvestment(2:end) HPFilteredRealWage(2:end) HPFilteredHours(2:end) Inflation ShortRate(2:end)]);
[T,N]=size(Y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[BLower,BUpper]=GetLimitsSmetsWouters(1);
[A,B]=GetPriorsSmetsWouters(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T0=500000;
rT=0.9;
Nt=5;
Ns=20;
MaxNEval=100000;
%
[Bhat,fopt,sa_nevals]=NKSimulAnnBayesianSmetsWouters(Y,T,BLower,BUpper,A,B,T0,rT,Nt,Ns,MaxNEval);
% For Euro area:
% load 'U:\CaiZhiqing\SeminarEmpiricalMacro\ModeLogPosterior'
% Bhat=[ ];
%
[YHAT,INST]=YhatSmetsWouters(Bhat,Y,T);
%
zz=1;
while zz<=N
    figure(1)
    subplot(2,N,zz)
    plot(Time,Y(:,zz),'k',Time,YHAT(:,zz),'r','LineWidth',2)
    xlim([Time(1) Time(end)])
    subplot(2,N,N+zz)
    plot(Time,zeros(size(Time)),'r',Time,Y(:,zz)-YHAT(:,zz),'k','LineWidth',2)
    xlim([Time(1) Time(end)])
    zz=zz+1;
end
% Here we get the Hessian:
%
P=GetPSmetsWouters(Bhat,Y,T);
[Bhat sqrt(diag(P))]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    Here we do Random Walk Metropolis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           I: The burn-in period
%
Desired=0.23; % Desired acceptance rate of the draws
BurnInNN=100000; % Number of burn-in iterations
SamplingInterval=1000;
ErgodicNN=1000*SamplingInterval; % Number of iterations of the ergodic distribution
IntervalForResettingOfc=10000;
%
c=0.5;
RR=c*chol(P)';
%
KK=size(Bhat,1);
BHAT=zeros(size(Bhat,1),BurnInNN);
DrawMean=zeros(size(Bhat,1),1);
DrawCovariance=eye(size(Bhat,1));
%
BB=Bhat;
% Log-likelihood conditional on the starting value:
LikeOld=-LogLikelihoodSmetsWouters(BB,Y,T);
BayesOld=LogPriorDensitySmetsWouters(BB,A,B,'Y');
%
% The number of accepted draws:
Accepted=0;
% The first element of the Markov chain:
BHAT(:,1)=BB;
%
jj=1;
while jj<BurnInNN
    % jj
    Index=-9999;
    while Index<0
        % Proposal draw:
        BBProposal=BB+RR*mvnrnd(DrawMean,DrawCovariance,1)';
        if sum(BBProposal>0)==KK
            LikeNew=-LogLikelihoodSmetsWouters(BBProposal,Y,T);
            if LikeNew~=-1.0e+100 & LikeNew~=LikeOld
                Index=1;
            end
        end
    end
    BayesNew=LogPriorDensitySmetsWouters(BBProposal,A,B,'Y');
    if BayesNew==0
        BHAT(:,jj)=BB;
        jj=jj+1;
    else
        R=exp(LikeNew-LikeOld)*(BayesNew/BayesOld);
        R=min([R 1]');
        r=rand;
        if R>=r
            Accepted=Accepted+1;
            BB=BBProposal;
            LikeOld=LikeNew;
            BayesOld=BayesNew;
        end
        BHAT(:,jj)=BB;
        jj=jj+1;
    end
    if (jj/IntervalForResettingOfc)==fix(jj/IntervalForResettingOfc)
        BurnInNN-jj
        Fraction=Accepted/IntervalForResettingOfc
        if Fraction<0.21
            RR=RR*0.9;
        elseif Fraction>0.25
            RR=RR*1.1;
        else
        end
        Accepted=0;
    end
end
% plot(BHAT(:,1:jj-1)')
% plot(BHAT(:,1:SamplingInterval:jj-1)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           II: The ergodic distribution:
clear DFILE varname BHAT
if Country=='US'
    DFILE='U:\CaiZhiqing\SeminarEmpiricalMacro\ErgodicDistributionUS';
elseif Country=='EA'
    DFILE='U:\CaiZhiqing\SeminarEmpiricalMacro\ErgodicDistributionEuroArea';
else
end
varname(1,:)='FRAC';
varname(2,:)='BHAT';
%
LikeOld=-LogLikelihoodSmetsWouters(BB,Y,T);
BayesOld=LogPriorDensitySmetsWouters(BB,A,B,'Y');
%
Accepted=0;
accepted=0;
Fraction=1;
%
jj=1;
while jj<ErgodicNN
    if (jj/10000)==fix(jj/10000)
        FRAC=Accepted/jj;
        ErgodicNN-jj
    end
    Index=-9999;
    while Index<0
        % Proposal draw:
        BBProposal=BB+RR*mvnrnd(DrawMean,DrawCovariance,1)';
        if sum(BBProposal>0)==KK
            LikeNew=-LogLikelihoodSmetsWouters(BBProposal,Y,T);
            if LikeNew~=-1.0e+100 & LikeNew~=LikeOld
                Index=1;
            end
        end
    end
    BayesNew=LogPriorDensitySmetsWouters(BBProposal,A,B,'Y');
    if BayesNew==0
        jj=jj+1;
    else
        R=exp(LikeNew-LikeOld)*(BayesNew/BayesOld);
        R=min([R 1]');
        r=rand;
        if R>=r
            accepted=accepted+1;
            Accepted=Accepted+1;
            BB=BBProposal;
            LikeOld=LikeNew;
            BayesOld=BayesNew;
        end
        jj=jj+1;
    end
    if (jj/SamplingInterval)==fix(jj/SamplingInterval)
        K=(jj/SamplingInterval)
        FRAC=Fraction;
        BHAT(:,K)=BB;
        save(DFILE,varname(1,:),varname(2,:));
    end
    if (jj/IntervalForResettingOfc)==fix(jj/IntervalForResettingOfc)
        Fraction=accepted/IntervalForResettingOfc
        if Fraction<0.21
            RR=RR*0.9;
        elseif Fraction>0.25
            RR=RR*1.1;
        else
        end
        accepted=0;
    end
end
FRAC=Accepted/ErgodicNN;
BHAT=BHAT';
save(DFILE,varname(1,:),varname(2,:));
%
[N,C]=size(BHAT);
SortedBHAT=sort(BHAT);
for jj=1:C
    Grid=linspace(min(SortedBHAT(:,jj)),max(SortedBHAT(:,jj)),100)';
    [PDF,Grid]=PDFContour(SortedBHAT(:,jj),Grid);
    [nn,Index]=max(PDF);
    Mode(jj,1)=Grid(Index);
end
[Mode SortedBHAT(fix(N*[0.05 0.95]'),:)']
MedianEstimate=SortedBHAT(fix(N*0.5),:)';
N=size(Y,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YHAT=YhatSmetsWouters(MedianEstimate,Y,T);
%
zz=1;
while zz<=N
    figure(1)
    subplot(1,N,zz)
    plot(Time,Y(:,zz),'k',Time,YHAT(:,zz),'r','LineWidth',1)
    xlim([Time(1) Time(end)])
    zz=zz+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Here we get the IRFs and fractions of FEV:
NN=size(BHAT,1);
Horizon=10*4;
%
Draw=1;
while Draw<NN
    [F,A0,H]=GetStateSpaceSmetsWouters(BHAT(Draw,:)',Y,T);
    FractionsOfFEV(:,:,:,Draw)=GetFEVsDSGEModel(F,A0,H,Horizon);
    IRFs(:,:,:,Draw)=GetIRFsDSGEModel(F,A0,H,Horizon);
    Draw=Draw+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Here we plot the IRFs and fractions of FEV:
HOR=(0:1:Horizon)';
%
Variable=1;
while Variable<=N
    Shock=1;
    while Shock<=N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Here we plot the IRFs:
        IRF=ExtractPercentiles(squeeze(IRFs(:,Variable,Shock,:))',[0.5 0.16 0.84 0.05 0.95]')';
        figure(2)
        subplot(N,N,(Shock-1)*N+Variable)
        plot(HOR,zeros(size(HOR)),'b:',HOR,IRF(:,1),'k',HOR,IRF(:,2:3),'r:',HOR,IRF(:,4:5),'r','LineWidth',2)
        xlim([0 Horizon])
        if Variable==1
            title('GDP')
        elseif Variable==2
            title('Consumption')
        elseif Variable==3
            title('Investment')
        elseif Variable==4
            title('Real wage')
        elseif Variable==5
            title('Hours')
        elseif Variable==6
            title('Inflation')
        elseif Variable==7
            title('Short rate')
        else
        end
        if Shock==1
            xlabel('e_g(t)')
        elseif Shock==2
            xlabel('e_b(t)')
        elseif Shock==3
            xlabel('e_i(t)')
        elseif Shock==4
            xlabel('e_a(t)')
        elseif Shock==5
            xlabel('e_p(t)')
        elseif Shock==6
            xlabel('e_w(t)')
        elseif Shock==7
            xlabel('e_r(t)')
        else
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Here we plot the fractions of FEV:
        FEV=ExtractPercentiles(squeeze(FractionsOfFEV(:,Variable,Shock,:))',[0.5 0.16 0.84 0.05 0.95]')';
        figure(3)
        subplot(N,N,(Shock-1)*N+Variable)
        plot(HOR,FEV(:,1),'k',HOR,FEV(:,2:3),'r:',HOR,FEV(:,4:5),'r','LineWidth',2)
        axis([0 Horizon 0 1])
        if Variable==1
            title('GDP')
        elseif Variable==2
            title('Consumption')
        elseif Variable==3
            title('Investment')
        elseif Variable==4
            title('Real wage')
        elseif Variable==5
            title('Hours')
        elseif Variable==6
            title('Inflation')
        elseif Variable==7
            title('Short rate')
        else
        end
        if Shock==1
            xlabel('e_g(t)')
        elseif Shock==2
            xlabel('e_b(t)')
        elseif Shock==3
            xlabel('e_i(t)')
        elseif Shock==4
            xlabel('e_a(t)')
        elseif Shock==5
            xlabel('e_p(t)')
        elseif Shock==6
            xlabel('e_w(t)')
        elseif Shock==7
            xlabel('e_r(t)')
        else
        end
        Shock=Shock+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    Variable=Variable+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

