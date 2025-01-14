function [Bhat,fopt,sa_nevals]=NKSimulAnnBayesianSmetsWouters(Y,T,BLower,BUpper,A,B,T0,rT,Nt,Ns,MaxNEval);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DFILE=['U:\CaiZhiqing\SeminarEmpiricalMacro\ModeLogPosterior']; % ==> Output file
DFILE=['D:\code_github\matlab\dsge'];
varname(1,:)=['Bhat'];
warning('off','all')
BLower=BLower(:);
BUpper=BUpper(:);
%
sa_neps=4; % number of times eps tolerance is achieved before termination
sa_eps=1e-6; % convergence criteria
sa_nargs=length(BLower); % number of parameters
sa_nobds=0;
sa_nacc=0; % number of acceptions
sa_nevals=0; % number of evaluations
sa_opteval=0; % optimum number of function evaluations
%
fstar=Inf*ones(sa_neps,1);
%
MinusLogLikelihood=1.0000e+100;
while MinusLogLikelihood==1.0000e+100
    x=BLower+(BUpper-BLower).*rand(sa_nargs,1); % Starting values for model parameters
    MinusLogLikelihood=LogLikelihoodSmetsWouters(x,Y,T);
end
f=MinusLogLikelihood+LogPriorDensitySmetsWouters(x,A,B);
sa_nevals=sa_nevals+1;
Bhat=x;
fopt=f;
xtot=x;
fstar(1)=f;
% Maximum step size:
VM=(BUpper-BLower);
% Loop starts:
while 1
    nup=0; % Number of uphill movements:
    nrej=0; % Number of rejections:
    nnew=0;
    ndown=0; % Number of downhill movements:
    lnobds=0;
    nacp=zeros(sa_nargs,1);
    for m=1:Nt
        for j=1:Ns
            % For each single parameter:
            for h=1:sa_nargs
                if sa_nevals>=MaxNEval
                    disp('Too many function evaluations')
                    return
                end
                % Generate xp, trial value of x:
                MinusLogLikelihood=1.0000e+100;
                while MinusLogLikelihood==1.0000e+100
                    xp=x;
                    xp(h)=x(h)+VM(h)*(2*rand(1,1)-1.0);
                    if (xp(h)<BLower(h)) | (xp(h)>BUpper(h))
                        xp(h)=BLower(h)+(BUpper(h)-BLower(h))*rand(1,1);
                        lnobds=lnobds+1;
                        sa_nobds=sa_nobds+1;
                    end
                    MinusLogLikelihood=LogLikelihoodSmetsWouters(xp,Y,T);
                end
                fp=MinusLogLikelihood+LogPriorDensitySmetsWouters(xp,A,B);
                if fp==-Inf
                    fp=Inf;
                end
                sa_nevals=sa_nevals+1;
                if fix(sa_nevals/1000)==(sa_nevals/1000)
                    save(DFILE,varname(1,:))
                    MaxNEval-sa_nevals
                end
                if fp<=f % Minimize: accept if the function value decreases:
                    x=xp;
                    f=fp;
                    sa_nacc=sa_nacc+1;
                    nacp(h)=nacp(h)+1;
                    nup=nup+1;
                    if fp<fopt % If smaller than any previous point, record as new optimum:
                        Bhat=xp;
                        fopt=fp;
                        sa_opteval=sa_nevals;
                        nnew=nnew+1;
                    end
                else % Function value increases:
                    p=exp((f-fp)/T0);
                    pp=rand(1,1); % Random number:
                    if pp<p
                        x=xp;
                        f=fp;
                        sa_nacc=sa_nacc+1;
                        nacp(h)=nacp(h)+1;
                        ndown=ndown+1;
                    else
                        nrej=nrej+1;
                    end
                end
            end
        end
        % Adjust maximal step size VM:
        c=ones(sa_nargs,1)*2;
        for i=1:sa_nargs
            ratio=nacp(i)/Ns;
            if ratio>0.6
                VM(i)=VM(i) * (1+c(i)*(ratio-0.6)/0.4);
            elseif ratio <0.4
                VM(i)=VM(i)/(1+c(i)*((0.4-ratio)/0.4));
            end
            if VM(i)>(BUpper(i)-BLower(i))
                VM(i)=BUpper(i)-BLower(i);
            end
        end
        for i=1:sa_nargs
            nacp(i) = 0;
        end
    end
    % Check termination criteria:
    fstar(1)=f;
    quit = ((fstar(1)-fopt) <= sa_eps);
    if any(abs(fstar-f)>sa_eps)
        quit=0;
    end
    if quit
        disp(['Simulated annealing achieved termination after ', num2str(sa_nevals),' evaluations']);
        return
    end
    % Reduce temperature:
    T0=T0*rT;
    fstar(2:4)=fstar(1:3);
    % Continue from current optimum:
    x=Bhat;
    f=fopt;
end
Bhat=Bhat(:);



