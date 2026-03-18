function MMIPDE
clear all;
close all;
clc;
fun = 't30';
foutput = 'MMIPDE';
prob='POE_4obj_nlinear';
foutput2 = strcat('Resultados_t30_',prob,'_MMIPDE');
nloop=100;
nsol=10;
nvar=4;
narchive=2*nsol;
nrun=10;
% For finding a Pareto front, xPareto ={{x1}, {x2},...}
% to minimise {f}={f1, f2, ...,fM}^T
% subject to {lb} <= {x} <= ub

% lower and upper boinds of {x}

lb = [0.1;0.1;0.1;0.1]; 
ub = [2;2;2;2];
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Wansasueb, K., Pholdee, N., Panagant, N., & Bureerat, S.
%   (2020). Multiobjective meta-heuristic with iterative parameter
%   distribution estimation for aeroelastic design of an aircraft wing.
%   Engineering with Computers, doi:10.1007/s00366-020-01077-w       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fun = objective function file name (string)
% fout = output data file name (string)
% nvar = number of design variables
% nloop = number of iterations
% nsol = population size
% narchive = external archive size
% xPareto = non-dominated soltions
% fPareto = objective function values of the solutions in xPareto
% gPareto = constraint values of the solutions in xPareto
% bPareto = non-dominated binary strings for {xh}
%%%%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%%%%

tic;
for k=1:nrun % No. of Run 
% Initialisation
rand('state',sum(100*clock));% reset randomisation (for older MATLAB versions)
nvar0=5;% no. of adaptive parameters {xh}
% For Fl, Fu, CRl, CRu, and MutationScheme 
a0=[1;0.75;0.1;0.75;1];% lower bounds of the parameters
b0=[2; 1.5;0.5; 1.0;5];% upper bounds of the parameters
nbit0=7;% no. of binary bits for estimating the parameters

[x0,f0,g0]=mmipde_initial(fun,nsol,nvar,lb,ub);
% x0
% randomly generated binary population for the control parameters
bin0 = round(rand(nbit0*nvar0,nsol));

% performe non-dominated sorting
[bPareto,xPareto,fPareto,gPareto,A]=mmipde_selection(bin0,x0,f0,g0,[],[],[],[],[],narchive);

nvector=ceil(nsol/10);% no. of probability vectors
Pmatrix=0.5*ones(nvector,nvar0*nbit0);% initial probaility matrix

nsol0=ceil(nsol/nvector)*ones(nvector);% this will be used for udating a probability matrix later
nsol0(nvector)=nsol-sum(nsol0(1:(nvector-1)));
  
   iter=1;% iteration counter
    while iter <= nloop
        clc;        
        fprintf('Run = %d\n',k); 
        fprintf('Ger = %d\n',iter); 
        iter=iter+1;
        % Generate a population of the control parameters 
        % according to the given
        % probability matrix Pmatrix and nsol
        bin0=[];
        for i=1:nvector
            bin0i=mmipde_findsubpop(Pmatrix(i,:),nsol0(i));
            bin0=[bin0 bin0i];
        end

        % convert binary to real control parameters {xh} in the paper
        xh=bin2real(bin0,a0,b0);

        % reproduction operator {xl} in the paper
        x1=mmipde_reproduct(x0,xPareto,lb,ub,xh);
        f1=[];g1=[];

        % function evaluations
        for i=1:nsol
%             x1 = ones(nvar,1);
            [f1(:,i),g1(:,i)]=feval(fun,x1(:,i));
        end

        % update the external Pareto archive
        [bPareto,xPareto,fPareto,gPareto,A]=mmipde_selection(bin0,x1,f1,g1,bPareto,xPareto,fPareto,gPareto,A,narchive);



        % provide bestsol used to update the probability matrix (Pmatrix)
        % Scheme1
        nac0=size(bPareto,2);
        nshff=randperm(nac0);
        bPareto0=bPareto(:,nshff);
        nbestsol=floor(nac0/nvector);
        nbs=[0 cumsum(nbestsol*ones(1,nvector))];
        if nbestsol < 1
            for i=1:nvector
                bestsol(:,i)=mean(bPareto0,2);
            end
        else
            for i=1:nvector
                bestsol(:,i)=mean(bPareto0(:,nbs(i)+1:nbs(i+1)),2);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%update the probability vector %%%%%%%%%
        Pmatrix0=Pmatrix;
        LR=0.5+0.1*rand*(-1)^round(rand);
        for i=1:nvector
            Pmatrix(i,:)=min(0.9,max(0.1,(1-LR)*Pmatrix(i,:) + LR*bestsol(:,i)'));
        end

    %     figure(2),clf
    %     subplot(3,1,1),hold on,bar(sort(round(xx0(1,:)))),xlabel('solution numbers'),ylabel('choices for F')
    %     plot([nsol/2 nsol/2],[0 2],'r-')
    %     subplot(3,1,2),bar(sort(round(xx0(5,:)))),xlabel('solution numbers'),ylabel('choices for mutation')
    %     subplot(3,1,3),bar(sort(xx0(4,:))),xlabel('solution numbers'),ylabel('CR')

        x0=x1;f0=f1;


    %     figure(1),clf,hold on
    %     plot(fPareto(1,:),fPareto(2,:),'d','markerfacecolor','g')
    %     plot(f1(1,:),f1(2,:),'or')
    %     title(['Iteration no. ' num2str(iter)])

    end
    final_candidatas = xPareto';
    final_objetivos = fPareto';
    for i=1:size(final_objetivos)
        aux(i,:)=[final_candidatas(i,:) final_objetivos(i,:)];
    end
    Execucao{k} = aux;
    resultados=strcat(foutput2);
    save(resultados,'Execucao');
end
runningTime=toc;
save(foutput,'bPareto','xPareto','fPareto','gPareto','runningTime')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%sub-programs%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x1 = mmipde_reproduct(x0,xPareto,a,b,xh)
% Algorithm 3 in the paper
[m0,n0]=size(x0);
[m1,n1]=size(xPareto);

npareto=size(xPareto,2);

for i=1:n0
    if npareto==1
        xbest1=xPareto;
        xbest2=xbest1;
    else
        nb=randperm(n1);
        xbest1=xPareto(:,nb(1));
        xbest2=xPareto(:,nb(2));
    end
    nr=randperm(n0);
    i1=nr(1);i2=nr(2);i3=nr(3);i4=nr(4);
    xr1=x0(:,i1);% Randomly seletced individual 1
    xr2=x0(:,i2);% Randomly seletced individual 2
    xr3=x0(:,i3);% Randomly seletced individual 3
    xr4=x0(:,i4);% Randomly seletced individual 4
    
    if round(xh(1,i))==1
        F=normrand(xh(2,i),xh(3,i));
    else
        F=unifrand(xh(3,i),xh(2,i));
    end
    Imutate=round(xh(5,i));CR=xh(4,i);
    switch Imutate
        case 1
            vi=EAreproduct(xbest1,xbest2,a,b,1-CR);
        case 2
            ui=x0(:,i)+F*(rand*xbest1+rand*xbest2-x0(:,i));
            vi=DEcrossover(ui,x0(:,i),CR,a,b);
        case 3
            ui=rand*xbest1+rand*xbest2+F*(rand*xr1...
                +rand*xr2-rand*xr3-rand*xr4);
            vi=DEcrossover(ui,x0(:,i),CR,a,b);
        case 4
            ui=x0(:,i)+F*(xbest1-x0(:,i));
            vi=DEcrossover(ui,x0(:,i),CR,a,b);
        case 5
            vi=Simplemutate(xbest1,a,b);
            %             ui=xr1+F*(xr2+(-1)^round(rand)*xr3);
            %             vi=DEcrossover(ui,x0(:,i),CR,a,b);
    end
    
    x1(:,i)=vi;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pr=normrand(mu0,sigma0)
pr=mu0+sigma0.*randn(size(mu0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pr=unifrand(lb,ub)
pr=lb+(ub-lb)*randn(size(lb));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v=DEcrossover(u,x,CR,a,b)
for i=1:length(u) % binomial crossover
    if rand < CR
        v(i,1)=max(a(i),min(b(i),u(i)));
    else
        v(i,1)=x(i);
    end
end
%%%%%%%%%%%%%%
function v=EAreproduct(x1,x2,a,b,pm)
n=length(x1);
stl=-0.25+1.5*rand;
for i=1:n
    PR=rand;
    if PR < 0.33
        if rand < 0.5
            v(i,1)=x1(i);
        else
            v(i,1)=x2(i);
        end
    elseif  PR >= 0.33 & PR < 0.66
        v(i,1)=x1(i)+stl*(x2(i)-x1(i));
    else
        v(i,1)=x1(i)+(-0.25+1.5*rand)*(x2(i)-x1(i));
    end
    
    if v(i,1) < a(i)
        v(i,1)=a(i);
    elseif v(i,1) > b(i)
        v(i,1)=b(i);
    end
end
if rand < pm
    v=Simplemutate(v,a,b);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v=Simplemutate(v,a,b)
n=length(v);
irow=max([ceil((rand)*n) 1]);
PR=rand;
if PR < 0.5
    v(irow,1)=min(b(irow),max(a(irow),a(irow)+1.15*rand*(v(irow,1)-a(irow))));
else
    v(irow,1)=min(b(irow),max(a(irow),v(irow,1)+1.15*rand*(b(irow)-v(irow,1))));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate nsol binary solutions from the given propability vector pop0
function bin0=mmipde_findsubpop(pop0,nsol)
[m,n]=size(pop0);
n1=round(nsol*pop0);
bin0=zeros(n,nsol);
for i=1:n;
    Ishuff=randperm(nsol);
    bin0(i,1:n1(i))=ones(1,n1(i));
    bin0(i,:)=bin0(i,Ishuff);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=bin2real(bin,a,b)
[m,n]=size(bin);
nvar=length(a);
nbit=m/nvar;

for i=1:n
    for j=1:nvar
        x(j,i)=bin2dec(bin((j-1)*nbit+1:j*nbit,i),a(j),b(j));
    end
end
function x=bin2dec(bin,a,b)

%
% Transformation from binary string to real number
% with lowr limit a and upper limit b

n=max(size(bin));
trans=cumprod(2*ones(size(bin)))/2;
real1=sum(bin.*trans);

x=a+(real1*(b-a))/(2^n-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update the external Pareto archive
function [paretob1,pareto1,fPareto1,gPareto1,A]=mmipde_selection(bin1,x1,f1,g1,paretob,pareto,fPareto,gPareto,A,narchive)
% figure(2),clf,hold on
% plot(f1(1,:),f1(2,:),'*b')
% plot(fPareto(1,:),fPareto(2,:),'*r')

bin=[paretob bin1];
x=[pareto x1];
f=[fPareto f1];
g=[gPareto g1];

[m0,n0]=size(fPareto);
[m1,n1]=size(x);

for i=1:n1
    fi=f(:,i);
    gi=g(:,i);
    A(i,i)=0;
    for j=(n0+1):n1
        fj=f(:,j);
        gj=g(:,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        [p_count1,p_count2]=fdominated(fi,gi,fj,gj);
        A(i,j)=p_count1;
        A(j,i)=p_count2;
        %%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=(n0+1):n1
    fi=f(:,i);
    gi=g(:,i);
    A(i,i)=0;
    for j=(i+1):n1
        fj=f(:,j);
        gj=g(:,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        [p_count1,p_count2]=fdominated(fi,gi,fj,gj);
        A(i,j)=p_count1;
        A(j,i)=p_count2;
        %%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=sum(A,1);
Indm=[];
for i=1:length(B)
    if B(i)==0
        Indm=[Indm i];
    end
end
nndm=length(Indm);
paretob1=bin(:,Indm);
pareto1=x(:,Indm);
fPareto1=f(:,Indm);
gPareto1=g(:,Indm);
A=A(Indm,Indm);

if nndm > narchive
    nsl=fTASI(fPareto1,narchive);
    paretob1=paretob1(:,nsl);
    pareto1=pareto1(:,nsl);
    fPareto1=fPareto1(:,nsl);
    gPareto1=gPareto1(:,nsl);
    A=A(nsl,nsl);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seletion used when the number of non-dominated solutions
% exceeds the archive size
% Truncation Algorithm with Similar Individuals from 
% Gao, J. and J. Wang, WBMOAIS: A novel artificial immune 
% system for multiobjective optimization. Computers & Operations
% Research, 2010. 37(1): p. 50-61. 
function iselect=fTASI(f,Nmem)
[nvar,Npop]=size(f);
Nremove=Npop-Nmem;
% Compute Eucledian distance
for i=1:Npop
    fi=f(:,i);
    D(i,i)=0;% Eucledian distance
    for j=(i+1):Npop;
        fj=f(:,j);
        D(i,j)=sqrt(sum((fi-fj).*(fi-fj)));D(j,i)=D(i,j);
    end
end
KK=1:Npop;
iremove=[];
while length(iremove)<Nremove
    % Find the closest pair
    Dmin=[];J1=[];
    for i=1:(size(D,2)-1)
        J1=(i+1):size(D,2);
        [Dmin(i),nmin]=min(D(i,J1));
        J(i)=J1(nmin);
    end
    
    [Dmin2,nmin2]=min(Dmin);
    I=nmin2;J=J(I);
    fI=f(:,I);fJ=f(:,J);
    DI=D(I,:);DI(I)=[];DI=sort(DI);
    DJ=D(J,:);DJ(J)=[];DJ=sort(DJ);
    sw=1;k=0;
    while sw==1
        k=k+1;
        if DI(k) < DJ(k)
            K=I;sw=0;
        elseif DJ(k) < DI(k)
            K=J;sw=0;
        end
        if k==length(DI)
            K=I;sw=0;
        end
    end
    iremove=[iremove KK(K)];
    f(:,K)=[];KK(K)=[];D(:,K)=[];D(K,:)=[];
end
iselect=setdiff(1:Npop,iremove);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constrained optimisation domination sorting
function [p1,p2]=fdominated(f1,g1,f2,g2)
n=length(f1);
mg1=max(g1);
mg2=max(g2);

icount11=0;
icount12=0;
icount21=0;
icount22=0;

if mg1<=0&mg2<=0
    for i=1:n
        if f1(i) <= f2(i)
            icount11=icount11+1;
        end
        if f1(i) < f2(i)
            icount12=icount12+1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        if f2(i) <= f1(i)
            icount21=icount21+1;
        end
        if f2(i) < f1(i)
            icount22=icount22+1;
        end
    end
    if icount11 == n & icount12 > 0
        p1=1;
    else
        p1=0;
    end
    if icount21 == n & icount22 > 0
        p2=1;
    else
        p2=0;
    end
elseif mg1 <=0 & mg2 > 0
    p1=1;p2=0;
elseif mg2 <=0 & mg1 > 0
    p1=0;p2=1;
else
    if mg1 <= mg2
        p1=1;p2=0;
    else
        p1=0;p2=1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,f,g] = mmipde_initial(fun,nsol,nvar,a,b)
%
% Randomly initiate the population, design variables
% nvar=no. of variables
% nbit is the number of cell in each variable
% nsol is a number of gene
%
for i=1:nvar
    x(i,:)=a(i)+(b(i)-a(i))*(randperm(nsol)-1+rand(1,nsol))/nsol;
end
for i=1:nsol
    [f(:,i),g(:,i)]=feval(fun,x(:,i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% End of file %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%