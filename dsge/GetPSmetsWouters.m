function P=GetPSmetsWouters(Bhat,Y,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=size(Bhat,1);
BhatMaxLogPosterior=Bhat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delta=1.0e-14;
DeltaSaved=Delta;
BhatSaved=Bhat;
Delta=abs(Bhat*Delta);
Index=find((Delta)<eps);
Delta(Index)=ones(size(Index))*eps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii=1;
CheckTHETA=0;
while ii<=N
    bhatL=Bhat;
    bhatL(ii)=bhatL(ii)-Delta(ii)/2;
    [~,TrendLLower]=LogLikelihoodSmetsWouters(bhatL,Y,T);
    if TrendLLower==-9999
        CheckTHETA=-9999;
        break
    end
    bhatU=Bhat;
    bhatU(ii)=bhatU(ii)+Delta(ii)/2;
    [~,TrendLLUpper]=LogLikelihoodSmetsWouters(bhatU,Y,T);
    if TrendLLUpper==-9999
        CheckTHETA=-9999;
        break
    end
    THETA(:,ii)=(TrendLLUpper-TrendLLower)/Delta(ii);
    ii=ii+1;
end
POP=zeros(size(Bhat,1));
for tt=1:size(THETA,1)
    theta=THETA(tt,:)';
    POP=POP+theta*theta';
end
P=MyInverse(POP);



