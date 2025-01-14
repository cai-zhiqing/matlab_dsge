function F=GetBetaParametersF(x,modeB,stdB)
Bhat=exp(x)+1;
alpha=Bhat(1);
beta=Bhat(2);
varB=stdB^2;
F=(modeB-(alpha-1)/(alpha+beta-2))^2+(varB-(alpha*beta)/((alpha+beta)^2*(alpha+beta+1)))^2;




