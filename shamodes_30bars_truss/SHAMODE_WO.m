% Copyright (c) 2019, Natee Panagant
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of  nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function rst = SHAMODE_WO(fun,nloop,nsol,nvar,narchive,a,b,run)
% This algorithm, Success History肪ased Adaptive Multi-Objective Differential Evolution
% with Whale Optimization (SHAMODE-WO) is an improved multiobjective version of Success
% History-based Adaptive Differential Evolution (SHADE) by integrating modified adaptive
% strategies, non-dominated sorting algorithm, and additional population update operator
% from Whale Optimization Algorithm (WOA).
%
% The algorithm is published in:
%
% Panagant, N., Bureerat, S., & Tai, K. (2019). A novel self-adaptive hybrid
% multi-objective meta-heuristic for reliability design of trusses with simultaneous
% topology, shape and sizing optimisation design variables. Structural and
% Multidisciplinary Optimization, 60(5), 1937-1955.
% DOI: https://doi.org/10.1007/s00158-019-02302-x
%
    rand('state',sum(100*clock));
    tic
    
    min_nsol=4;%
    max_nsol=nsol;
    current_nsol=nsol;


    % Generate initial population
    [x0,f0,g0] = moea_initialisation(fun,nsol,nvar,a,b)
    A=zeros(nsol);
    A0=A;
    [ppareto,fpareto,gpareto,A]=pareto_sorting(x0,f0,g0,[],[],[],A,nsol);
    
    %Adaptive Processor (Initialize)
    adapt_dat=struct();
    adapt_dat.achieve_ratio=1.4;
    adapt_dat.achieve_size=adapt_dat.achieve_ratio*current_nsol;
    adapt_dat.x0a=[];
    adapt_dat.f0a=[];
    adapt_dat.g0a=[];
%     adapt_dat.xpbest_ratio=0.11;
    adapt_dat.mem_size=5;
    adapt_dat.mem_pos=1;
    adapt_dat.mem_F=0.5*ones(1,adapt_dat.mem_size);
    adapt_dat.mem_CR=0.5*ones(1,adapt_dat.mem_size);
    
    % MAIN LOOP
    for iter=1:nloop
        clc;        
        fprintf('Run = %f\n',run); % Nome do Algoritmo
        fprintf('Ger = %f\n',iter); % Nome do Algoritmo
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    DE operatos & FSD resizing   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update population
        [x1,f1,g1,adapt_dat]=de_reproduct(fun,x0,f0,g0,ppareto,fpareto,gpareto,a,b,iter,nloop,adapt_dat);
        % Selection
        [x2,f2,g2,A2,nsort]=DEselection(x1,f1,g1,x0,f0,g0,A0);
        [ppareto,fpareto,gpareto,A]=pareto_sorting(x2,f2,g2,ppareto,fpareto,gpareto,A,narchive);

        x0=x2;f0=f2;g0=g2;A0=A2;
    
        %DE Adaptive Processor
        [adapt_dat]=DE_adaptive_processor(adapt_dat,x0,f0,g0,nsort,nsol);

        % Save Results
        rst.ppareto{iter}=ppareto;
        rst.fpareto{iter}=fpareto;
        rst.gpareto{iter}=gpareto;
        rst.timestamp=datetime('now');
    end
end

%%%%%%%%%%%%%%%
% SUB-Program %
%%%%%%%%%%%%%%%
function [x1,f1,g1,adapt_dat] = de_reproduct(fun,x0,f0,g0,ppareto,fpareto,gpareto,a,b,iter,nloop,adapt_dat)
    [nvar,nsol]=size(x0);
    f1=zeros(size(f0,1),nsol);
    g1=zeros(size(g0,1),nsol);
    
    % DE point generation for n1 parents (x1)
    % x0 = population from the previous generation
    % f0 = corresponding objective values
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % DE Mutation Operator %
    %%%%%%%%%%%%%%%%%%%%%%%%
    %Generate scaling factror(F)
    F=zeros(1,nsol);
    F_ind=ceil(adapt_dat.mem_size*rand(1,nsol));
    F_mid=adapt_dat.mem_F(F_ind);
    pos=true(size(F));
    while sum(pos)>0
        % Random Scaling Factor from Cauchy Distribution
        F(pos)=F_mid(pos)+0.1*tan(pi*(rand(1,sum(pos))-0.5)); % Eq. (3) in the published paper
        pos=F<=0;
    end
    F=min(1,F);
    adapt_dat.F=F;
    x0a=[x0,adapt_dat.x0a];
    
    Npbest=size(ppareto,2);
    pbi=ceil(Npbest*rand(1,nsol));
    xpbest=ppareto(:,pbi);
    
    [i1,i2]=ind_gen(nsol,size(x0a,2));
    xm=x0+F(ones(nvar,1),:).*(xpbest-x0+x0(:,i1)-x0a(:,i2));   % Eq. (1) in the published paper

    xm=LUbound(xm,x0,[a,b]);
    % Updatae with Bubble-net attacking method
    % of Whale Optimization Algorithm (WOA)
    bb=1;
    pbi2=ceil(Npbest*rand(1,nsol));
    xpbest2=ppareto(:,pbi2);
    for wi=1:nsol
        if rand<0.5
            % Eq. (7-8) in the published paper
            l=2*rand-1;
            xm(:,wi)=abs(xpbest2(:,wi)-xm(:,wi))*exp(bb*l)*cos(2*pi*l)+xpbest2(:,wi);
        end
    end
    xm=LUbound(xm,x0,[a,b]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % DE Crossover Operation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Generate crossover ratio
    CR_ind=ceil(adapt_dat.mem_size*rand(1,nsol));
    CR_mid=adapt_dat.mem_CR(CR_ind);
    CR=normrnd(CR_mid, 0.1);    % Eq. (4) in the published paper
    CR=max(0,min(1,CR));
    adapt_dat.CR=CR;
    
    % Binomial Crossover (Line 115-121 is Eq. (2) in the published paper)
    x1=x0;
    mask=rand(nvar,nsol)<=CR(ones(nvar,1),:);%mask is used to indicate which elements of xr comes from xm
    rows=floor(rand(1,nsol)*nvar)+1; % choose one position where the element of xr always come from xm
    cols=(1:nsol);
    jrand=sub2ind([nvar,nsol],rows,cols);
    mask(jrand)=true;
    x1(mask)=xm(mask);
    %
    for i=1:size(x1,2)
        [f1(:,i),g1(:,i)]=feval(fun,x1(:,i));
    end
end

function [adapt_dat] = DE_adaptive_processor(adapt_dat,x0,f0,g0,nsort,nsol)
    % Get success update index (Parent of success offspring)
    % index 1:nsol = Parent, (nsol+1):(2*nsol) = Offspring
    sind=nsort(1:nsol); % Success Offspring = Offspring in first 1:nsol rank
    sind=sind(sind>nsol)-nsol; % Offspring index - nsol = Parent index
    %
    adapt_dat.achieve_size=round(adapt_dat.achieve_ratio*nsol);

    % Update adaptive population
    x0a=[adapt_dat.x0a,x0(:,sind)];
    f0a=[adapt_dat.f0a,f0(:,sind)];
    g0a=[adapt_dat.g0a,g0(:,sind)];
    [~,Uind]=unique(x0a','rows');
    x0a=x0a(:,Uind);
    f0a=f0a(:,Uind);
    g0a=g0a(:,Uind);
    
    if size(adapt_dat.f0a,2)<=adapt_dat.achieve_size
        adapt_dat.x0a=x0a;
        adapt_dat.f0a=f0a;
        adapt_dat.g0a=g0a;
    else
        ind=randperm(size(adapt_dat.f0a,2));
        ind=ind(1:adapt_dat.achieve_size);
        adapt_dat.x0a=x0a(:,ind);
        adapt_dat.f0a=f0a(:,ind);
        adapt_dat.g0a=g0a(:,ind);
    end

    %Update scaling factror and crossover memory
    if sum(sind)>0 && adapt_dat.mem_CR(adapt_dat.mem_pos)~=-1
        %Update F
        sF=adapt_dat.F(sind);
        lmeanF=sum(sF.^2)/sum(sF); % Lehmer mean Eq. (5) in the published paper
        adapt_dat.mem_F(adapt_dat.mem_pos,:)=lmeanF;    

        %Update CR
        sCR=adapt_dat.CR(sind);
        if max(sCR)==0 || adapt_dat.mem_CR(adapt_dat.mem_pos)==-1
            adapt_dat.mem_CR(adapt_dat.mem_pos)=-1;
        else
            lmeanCR=sum(sCR.^2)/sum(sCR); % Lehmer mean Eq. (6) in the published paper
            adapt_dat.mem_CR(adapt_dat.mem_pos)=lmeanCR;
        end
    end
    adapt_dat.mem_F=max(0,min(1,adapt_dat.mem_F));
    adapt_dat.mem_CR=max(0,min(1,adapt_dat.mem_CR));
    
    %Shift memory position
    adapt_dat.mem_pos=adapt_dat.mem_pos+1;
    if adapt_dat.mem_pos>adapt_dat.mem_size
        adapt_dat.mem_pos=1;
    end
end

function [i1,i2] = ind_gen(nsol,nsolA)
    % Generate no duplicate index for mutation process
    i0=1:nsol;

    i1=max(1,ceil(rand(1,nsol)*nsol));
    dpos=i1==i0;
    iter=0;
    while sum(dpos)>0
        iter=iter+1;
        i1(dpos)=max(1,ceil(rand(1,sum(dpos))*nsol));
        dpos=i1==i0;
        if iter>999
            error('check');
        end
    end
    
    i2=max(1,ceil(rand(1,nsol)*nsolA));
    dpos=i2==i0 | i2==i1;
    iter=0;
    while sum(dpos)>0
        iter=iter+1;
        i2(dpos)=max(1,ceil(rand(1,sum(dpos))*nsolA));
        dpos=i2==i0 | i2==i1;
        if iter>999
            error('check');
        end
    end
end

function [vi] = LUbound(vi,pop,lu)
    bound_ratio=2;
    lb=lu(:,ones(1,size(vi,2)));
    ub=lu(:,2*ones(1,size(vi,2)));
    lind=vi<lb;
    vi(lind)=(1/bound_ratio)*((bound_ratio-1)*lb(lind)+pop(lind));
    uind=vi>ub;
    vi(uind)=(1/bound_ratio)*((bound_ratio-1)*ub(uind)+pop(uind));
end

function [pop3,f3,g3,A3,nsort]=DEselection(pop2,f2,g2,pop1,f1,g1,A1)
    % DE selection for multiobjective optimization problems
    % First NP solutions with lowest level of being dominated are survived
    % NP = Number of population
    %
    pop=[pop1 pop2];
    f=[f1 f2];
    g=[g1 g2];
    
    [m0,n0]=size(A1);
    [m1,n1]=size(pop);
    [m2,n2]=size(f);
    [m3,n3]=size(f2);
    %%%%%%%%%%%%%%%%%%%%%%
    A=zeros(size(pop,2));
    for i=1:n0
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
        end
    end
    
    B=ndlevel_sort(A);% level of being dominated, 1 = non-dominated
    
    [B,nsort]=sort(B);

    pop=pop(:,nsort);
    f=f(:,nsort);
    g=g(:,nsort);
    A=A(nsort,nsort);

    nlevel=zeros(1,max(B));
    for i=1:n1;nlevel(B(i))=nlevel(B(i))+1;end
    nlevel2=cumsum(nlevel);

    ncheck=nlevel(1);icheck=1;
    while ncheck < n3
        icheck=icheck+1;
        ncheck=ncheck+nlevel(icheck);
    end

    if ncheck == n3
        pop3=pop(:,1:n3);
        f3=f(:,1:n3);
        g3=g(:,1:n3);
        A3=A(1:n3,1:n3);
    elseif ncheck > n3 && icheck==1
        nc2=1:nlevel2(1);

        popc=pop(:,nc2);
        fc=f(:,nc2);
        gc=g(:,nc2);
        Ac=A(nc2,nc2);

        iselect=farchive(fc,n3);

        pop3=popc(:,iselect);
        f3=fc(:,iselect);
        g3=gc(:,iselect);
        A3=Ac(iselect,iselect);
    else
        nn1=nlevel2(icheck-1);
        nn2=n3-nn1;

        nc1=1:nn1;
        nc2=(nn1+1):nlevel2(icheck);

        popc=pop(:,nc2);
        fc=f(:,nc2);
        gc=g(:,nc2);
        Ac=A(nc2,nc2);

        iselect=farchive(fc,nn2);

        pop3=pop(:,iselect);
        f3=f(:,iselect);
        g3=g(:,iselect);
        A3=A(iselect,iselect);

        nc3=[nc1,nn1+iselect];

        pop3=pop(:,nc3);
        f3=f(:,nc3);
        g3=g(:,nc3);
        A3=A(nc3,nc3);
    end
end
function [pareto1,fpareto1,gpareto1,A]=pareto_sorting(x1,f1,g1,pareto,fpareto,gpareto,A,narchive)
    % Merge Current population and Pareto front
    % Re-sort the merged solutions to generate new Pareto front
    %
    x=[pareto x1];
    f=[fpareto f1];
    g=[gpareto g1];

    [m0,n0]=size(fpareto);
    [m1,n1]=size(x);

    for i=1:n1
        xi=x(:,i);
        fi=f(:,i);
        gi=g(:,i);
        A(i,i)=0;
        for j=(n0+1):n1
            xj=x(:,j);
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
        xi=x(:,i);
        fi=f(:,i);
        gi=g(:,i);
        A(i,i)=0;
        for j=(i+1):n1
            xj=x(:,j);
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
    for i=1:n1
        if B(i)==0
            Indm=[Indm i];
        end
    end
    nndm=length(Indm);

    pareto1=x(:,Indm);
    fpareto1=f(:,Indm);
    gpareto1=g(:,Indm);
    A=A(Indm,Indm);

    if nndm > narchive
        nsl=farchive(fpareto1,narchive);
        pareto1=pareto1(:,nsl);
        fpareto1=fpareto1(:,nsl);
        gpareto1=gpareto1(:,nsl);
        A=A(nsl,nsl);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seletion used when the number of non-dominated solutions
% exceeds the archive size
function iselect=farchive(f,narchive)
    iselect=farchive22(f,narchive);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B=ndlevel_sort(A)
    [m,n]=size(A);
    b=sum(A,1);
    [b,ns]=sort(b);

    ilevel=1;
    B(1)=ilevel;

    for i=2:n
        if b(i)==b(i-1)
            B(i)=ilevel;
        else
            ilevel=ilevel+1;
            B(i)=ilevel;
        end
    end
    B(ns)=B;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p1,p2]=fdominated(f1,g1,f2,g2)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Non-dominated sorting %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    n=length(f1);
    mg1=max(g1);
    mg2=max(g2);

    icount11=0;
    icount12=0;
    icount21=0;
    icount22=0;

    if mg1<=0 && mg2<=0
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
        if icount11 == n && icount12 > 0
            p1=1;
        else
            p1=0;
        end
        if icount21 == n && icount22 > 0
            p2=1;
        else
            p2=0;
        end
    elseif mg1 <=0 && mg2 > 0
        p1=1;p2=0;
    elseif mg2 <=0 && mg1 > 0
        p1=0;p2=1;
    else
        if mg1 <= mg2
            p1=1;p2=0;
        else
            p1=0;p2=1;
        end
    end
end

function [x,f,g] = moea_initialisation(fun,nsol,nvar,a,b)
    % Generate initial solutions
    x=a+rand(nvar,nsol).*b;
    for i=1:nsol
        [f(:,i),g(:,i)]=feval(fun,x(:,i));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iselect=farchive22(f,narchive)
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Clustering Technique %
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Remove exceed solutions when number of solutions in Pareto front is
    % exceed archive size with clustering technique
    %
    [m,n]=size(f);

    fmax=max(f,[],2);fmin=min(f,[],2);
    fdel=max(fmax-fmin,1e-5);
    for i=1:n
        f(:,i)=(f(:,i)-fmin)./fdel;
    end
    nn=narchive;
    x1=linspace(0,1,nn);
    x2=1-x1;
    iselect0=1:n;
    for i=1:narchive
        df=[f(1,:)-x1(i);f(2,:)-x2(i)];
        d2=std(df);
        [d2min,imin]=min(d2);
        iselect(i)=iselect0(imin);
        iselect0(imin)=[];
        f(:,imin)=[];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%