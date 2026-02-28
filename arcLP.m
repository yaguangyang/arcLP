%
%                       arcLP.m Version 1.0
%
%                          Yaguang Yang
%
%
%                          August 2016
%
%
% arcLP.m implements an arc-search infeasible interior-point 
% algorithm for standard Linear Programming problems proposed in 
% [Y. Yang] Two computationally efficient polynomial-iteration 
% infeasible interior-point algorithms for linear programming
% Numerical Algorithms, (2018) 79:957–992. According to the numerical 
% test reported in [Y. Yang], this algorithm is competitive to 
% the popular Mehrotra's algorithm. Unlike Mehrotra's algorithm,
% arcLP searches the optimizer along an arc. The implemented algorithm
% is proved to converge with polynomial bound of O(\sqrt{n}L) 
% (see [Y. Yang]). Bug report should be sent to the author. 
%
% Synopsis:
%            function [x,obj,kk,infe]=arcLP(A,b,c,d) 
% Input arguments:
%   A: input matrix (real) with dimension nxn.
%   b: input vector (real) with dimension nx1.
%   c: input vector (real) with dimension nx1.
%   d: d=1 use a function to reduce A to full rank. d=1 is suggested
%          for degenerate problems to make matrix A full rank 
%      d=0 or without d will not reduce A to full rank
%   
%
% Output arguments:
%   x:    the optimal solution of the linear programming problem.
%   obj:  the optimum of the linear programming problem.
%   kk:   total iteration numbers used to find the optimal solution.
%   infe: the norm of the error || Ax-b ||.
%
%                   Conditions for External Use
%                   ===========================
%
%   1. Due  acknowledgment  must  be  made of the  use of this code in 
%      research  reports  or  publications.  Whenever such reports are 
%      released for  public  access, a copy should be forwarded to the 
%      author.
%   2. This code may only be used for research and development, unless 
%      it has been agreed  otherwise with the author in writing.

function [x,obj,kk,infe]=arcLP(A,b,c,d) 

if exist('d','var') == 0
    d = 0;
end
epsilon=0.00000001; 
rho=0.1;
sigma_min=0.000001; sigma_max=0.3;
maxAlpha2=0.99*pi/2;      % scaling factor 2
c_orig=c; %A_orig=A; b_orig=b; 
x_indx=1:1:length(c);
x_final=zeros(length(c),1);   % store final x_opt
lxfinal=length(x_final);
A_post=sparse([]);               % store info to recover x eliminated in rule 9
A_pst=sparse(zeros(1,lxfinal));  % one row of A_post
b_post=[];                       % store info to recover x eliminated in rule 9
x_tobe=[];                       % remember the index of x to be stored
%----------------------------
%check the problem conditions
%----------------------------
[m,n]=size(A);
totNonzero=nnz(A);

[minApre,maxApre]=calRatioCondition(A);
Aratio=maxApre/minApre;
%--------------
%preprocessing
%--------------
prepro=1;
csum=0;

while prepro==1
    prepro=0;
    
    %rule 1: romove empty row
    for i=m:-1:1 
        nz=nnz(A(i,:));
        if nz==0 
            if abs(b(i))<=epsilon    % if no error b(i)=0
                A(i,:)=[]; b(i)=[]; m=m-1; %lambda(i)=[];
                prepro=1;
            else
                error('infeasible')
            end
        end
    end
    disp('Rule 1: remove empty rows')
    disp([num2str(m),'  ', num2str(n)])

    %rule 3: remove empty column
    for i=n:-1:1
        nz=nnz(A(:,i));
        if nz==0
            if c(i)>=-epsilon*100 %if there is no error c(i)>=0
                A(:,i)=[]; c(i)=[]; x_indx(i)=[]; %s(i)=[]; x(i)=[];
                n=n-1;
                prepro=1;
            else
                error('infinite xi')
            end
        end
    end
    disp('Rule 3: remove empty columns')
    disp([num2str(m),'  ', num2str(n)])

    %rule 5: remove row of singleton. Identify problem "Scorpion" as infeasible
    %on my old computer
    i=m;
    while i>=1
        if nnz(A(i,:))==1 %find a singleton
            idj=find(A(i,:));
            xj=b(i)/A(i,idj); x_final(x_indx(idj))=xj;
            if xj<0
                error('the problem is not feasible')
            end
            csum=csum+c(idj)*xj;
            for k=1:length(b)
                b(k)=b(k)-A(k,idj)*xj;
            end
            A(:,idj)=[]; c(idj)=[]; x_indx(idj)=[]; %s(j)=[]; x(j)=[];
            A(i,:)=[]; b(i)=[]; %lambda(i)=[];
            n=n-1; m=m-1; prepro=1; %i=i-1;
        end
        i=i-1;
    end 
    disp('Rule 5: remove rows of singleton')
    disp([num2str(m),'  ', num2str(n)])

    % rule 7: remove fixed variables
    % b(i)=0 but all A(i,:)>=0 or A(i,:)<=0
    for i=m:-1:1
        if b(i)==0
            if nnz(A(i,:))==0 %min(A(i,:))==0 & max(A(i,:))==0
                A(i,:)=[]; b(i)=[];
                m=m-1; prepro=1;
            elseif min(A(i,:))==0 && max(A(i,:))>0
                idxa=find(A(i,:)>0);
                for j=length(idxa):-1:1
                    A(:,idxa(j))=[]; c(idxa(j))=[]; n=n-1; x_indx(idxa(j))=[];
                end
                A(i,:)=[]; b(i)=[]; %lambda(i)=[];
                m=m-1; prepro=1;
            elseif min(A(i,:))<0 && max(A(i,:))==0
                idxa=find(A(i,:)<0);
                for j=length(idxa):-1:1
                    A(:,idxa(j))=[]; c(idxa(j))=[]; n=n-1; x_indx(idxa(j))=[];
                end
                A(i,:)=[]; b(i)=[]; %lambda(i)=[];
                m=m-1; prepro=1;
            end
        elseif b(i)<0
            if min(A(i,:))==0
%                 [m n]=size(A)
                error('infeasible problem')
            end
        elseif b(i)>0
            if max(A(i,:))==0
%                 [m n]=size(A)
                error('infeasible problem')
            end
        end
    end
    disp('Rule 7: remove fixed variables')
    disp([num2str(m),'  ', num2str(n)])
    
      %rule 9: for row i, one A(i,k) has different sign from all other A(i,j)
      %and A(i,k) has the same sign as b(i)
      %A_post(:,n) and b_post exist; A_pst=sparse(1,n) and x_tobe also exist
    if Aratio<10^6  % why did I put this?
         i=m;
         while i>=1
             if b(i)>0
                 idxp=find(A(i,:)>0);
             elseif b(i)<0
                 idxp=find(A(i,:)<0);
             else
                 idxp=find(A(i,:)>0);
                 idxm=find(A(i,:)<0);
                 if length(idxm)<length(idxp)
                     tmp=idxp; idxp=idxm; idxm=tmp; clear tmp
                 end
             end
             if length(idxp)==1 
                 aip=A(i,idxp);
                 A_pst(x_indx(idxp))=1; b_post=[b_post;b(i)/aip]; 
                 % only one element in this row is positive
                 x_tobe=[x_tobe;x_indx(idxp)];
                 if b(i)>0 || (b(i)==0 && aip>0)
                     idxm=find(A(i,:)<0);
                 elseif b(i)<0 || (b(i)==0 && aip<0)
                     idxm=find(A(i,:)>0);
                 end
                 for jj=1:length(idxm)
                     A_pst(x_indx(idxm(jj)))=A(i,idxm(jj))/aip;
                 end
                 A_post=[A_post;A_pst];
                 A_pst=zeros(1,lxfinal);
                 for j=m:-1:1
                     if j~=i
                         b(j)=b(j)-A(j,idxp)*b(i)/aip;
                         for k=1:length(idxm)
                             A(j,idxm(k))=A(j,idxm(k))-A(j,idxp)*A(i,idxm(k))/aip;
                         end
                     end
                 end
                 csum=csum+c(idxp)*b(i)/aip;
                 for k=1:length(idxm)
                     c(idxm(k))=c(idxm(k))-c(idxp)*A(i,idxm(k))/aip;
                 end
                 c(idxp)=[]; A(i,:)=[]; A(:,idxp)=[]; b(i)=[]; x_indx(idxp)=[];
                 m=m-1; n=n-1; prepro=1; %i=i-1; 
             end
             i=i-1;
         end
         disp('Rule 9: solving positive variabls')
         disp([num2str(m),'  ', num2str(n)])
    end
end
[m1,n1]=size(A);

%-----------------
% make A full rank
%-----------------
if d==1
    fullRank=0;
    A0=A; b0=b; c0=c;
    [A,b,c]=makeAfull(A,b,c);
    [m,n]=size(A);
    if m1==m
        A=A0; b=b0; c=c0; 
        clear A0 b0 c0
        disp('original A is full rank')
    else
        fullRank=0;
        disp('original A is not full rank')
    end
end

%-------------------------------------------
%initial point (Lustig's method is the best)
%-------------------------------------------
e=ones(n,1); AAT=A*A'; % invAA=inv(A*A');
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merhodra's initial point
%%%%%%%%%%%%%%%%%%%%%%%%%%

ty=AAT\(A*c);  %ty=invAA*(A*c);
ts=c-A'*ty;
tx=A'*(AAT\b);  %tx=A'*(invAA*b);
dx=max(-1.1*min(tx),0);
ds=max(-1.1*min(ts),0);
sumx=0; sums=0;
for i=1:n
    sumx=tx(i)+dx;
    sums=ts(i)+ds;
end
tdx=dx+0.5*(tx+dx*e)'*(ts+ds*e)/sumx;
tds=ds+0.5*(tx+dx*e)'*(ts+ds*e)/sums;
xx=tx+tdx*e; 
ll=ty; 
ss=ts+tds*e; 
max1=max(max(norm(A*xx-b),norm(A'*ll+ss-c)),xx'*ss/n);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lustig's initial point
%%%%%%%%%%%%%%%%%%%%%%%%%%
norm1b=max(norm(b,1)/100, 100); tx=A'*(AAT\b);  %tx=A'*(invAA*b); 
tmpx1=max(-min(tx),norm1b); tmpy1=1+norm(c,1);
lambda=zeros(m,1);  x=zeros(n,1); s=zeros(n,1); %lambdaIdx=zeros(m,1);
for i=1:n
    x(i)=max(tmpx1, tx(i));
    if c(i)>tmpy1
        s(i)=c(i)+tmpy1;
    elseif c(i)<-tmpy1
        s(i)=-c(i);
    elseif c(i)>=0 && c(i)<tmpy1
        s(i)=c(i)+tmpy1;
    elseif -tmpy1<=c(i) && c(i)<=0
        s(i)=tmpy1;
    end
end

max2=max(max(norm(A*x-b),norm(A'*lambda+s-c)),x'*s/n);
if max1<max2
    x=xx; s=ss; lambda=ll;
    disp('Merhodra initial point is used')
else 
    disp('Lustig initial point is used')
end


mu=s'*x/n; rB=A*x-b; rC=A'*lambda+s-c;
kk=0; %iteration count
singleEarly=0; nu_k=1; 
sOld=s; xOld=x; lOld=lambda;
while (norm(rB)/max(1, norm(b)))+(norm(rC)/max(1, norm(c))) ...
        +(abs(c'*x-b'*lambda)/max(norm(b),max(1,norm(c))))>epsilon
    if d==1
        if fullRank==0 || singleEarly==1
            disp('make A full rank')
            xid=find(x<10^(-6));
            % when x->0, remove corresponding columns in A, c, x, s, and rC
            % if m<n, make A full matrix using Gaussian elimination
%             if length(xid)>0 
            if ~isempty(xid)
                for i=length(xid):-1:1
                    A(:,xid(i))=[]; c(xid(i))=[]; x_indx(xid(i))=[];
                    x(xid(i))=[]; s(xid(i))=[];
                    rC(xid(i))=[];
                    xOld(xid(i))=[]; sOld(xid(i))=[];
                end
                [A,b,c,x,s,lambda,rB,rC,lOld,xOld,sOld] = ...
                    makeAfull(A,b,c,x,s,lambda,rB,rC,lOld,xOld,sOld);
                [m,n]=size(A);
                e=ones(n,1);
%                 invS=sparse(diag(1./s)); 
%                 XinvS=sparse(diag(x.*(1./s))); 
%                 AXinvS=sparse(A*XinvS);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Algo 5.1 Step 1 
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    %------------------------------------------------------------------
    % dx ds dlambda
    %------------------------------------------------------------------
    invS=sparse(diag(1./s)); 
    XinvS=sparse(diag(x.*(1./s))); 
    AXinvS=sparse(A*XinvS);
    tmpA=sparse(AXinvS*A');
    tmpb=AXinvS*rC-b;
    dl=tmpA\tmpb;
    ds=rC-A'*dl;
    dx=x-XinvS*ds;

    %----------------------------
    %pl,ps,px,ql,qs,qx 
    %----------------------------
    tmpb=-A*invS*e*mu;
    pl=tmpA\tmpb;
    ps=-A'*pl;
    px=invS*e*mu-XinvS*ps;
    tmpb=2*A*invS*(dx.*ds);
    ql=tmpA\tmpb;
    qs=-A'*ql;
    qx=-XinvS*qs-2*invS*(dx.*ds); 
    clear tmpA tmpb Abackup bbackup
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Algo 5.1 Step 2 (x,s)(alpha,simga)>0
    %%%%%%%%%%%%%%%%%%%%%%%%

    %----------------------
    % Algo 6.1 select alpha & sigma 
    %----------------------
    sigma_lb=sigma_min; sigma_ub=sigma_max;
    while sigma_ub-sigma_lb>sqrt(epsilon)
        %----------------
        % Algo 6.1 Step 1 
        %----------------
        Sigma=sigma_lb+0.5*(sigma_ub-sigma_lb);
        %----------------
        % Algo 6.1 Step 2 
        %----------------
        xpp=Sigma*px+qx; xp=dx;
        spp=Sigma*ps+qs; sp=ds;
        alphaX=pi/2; alphaS=pi/2; alpha=pi/2;
        alphaXi=pi/2; alphaSi=pi/2;
        alphaMinus=pi/2; %Initial min_{px<0,ps<0} { min(alphaXi,alphaSi) }
        alphaPlus=pi/2;  %Initial min_{px>0,ps>0} { min(alphaXi,alphaSi) }
        for i=1:n
            phi_i=min(rho*min(x),nu_k);
            %-------------------
            % x(sigma,alphaS)>=phi_i
            %-------------------
            b1=x(i)-phi_i;
            b2=xp(i); alphaXi=pi/2;
            if b2==0 && b1+xpp(i) < 0 % Case 1
                alphaXi=acos((b1+xpp(i))/xpp(i));
                if alphaXi < alphaX 
                    alphaX=alphaXi;
                    if alphaX < 0
                        warning('check');
                    end
                end
            end
            if xpp(i)==0 && b2>b1  % Case 2
                alphaXi=asin((b1)/b2);
                if alphaXi < alphaX 
                    alphaX=alphaXi;
                    if alphaX < 0
                        warning('check');
                    end
                end
            end
            if b2>0 && xpp(i)>0 && b1+xpp(i) < sqrt(b2^2+xpp(i)^2) %Case 3
                alphaXi=asin((b1+xpp(i))/sqrt(b2^2+xpp(i)^2)) ...
                    -asin((xpp(i))/sqrt(b2^2+xpp(i)^2));
                if alphaXi < alphaX 
                    alphaX=alphaXi;
                    if alphaX < 0
                        warning('check');
                    end
                end
            end
            if b2>0 && xpp(i)<0 && b1+xpp(i) < sqrt(b2^2+xpp(i)^2) %Case 4
                alphaXi=asin((b1+xpp(i))/sqrt(b2^2+xpp(i)^2)) ...
                    +asin((-xpp(i))/sqrt(b2^2+xpp(i)^2));
                if alphaXi < alphaX 
                    alphaX=alphaXi;
                    if alphaX < 0
                        warning('check');
                    end
                end
            end
            if b2<0 && xpp(i)<0 && b1+xpp(i) < 0 %Case 5
                alphaXi=pi-asin(-(b1+xpp(i))/sqrt(b2^2+xpp(i)^2)) ...
                    -asin((-xpp(i))/sqrt(b2^2+xpp(i)^2));
                if alphaXi < alphaX 
                    alphaX=alphaXi;
                    if alphaX < 0
                        warning('check');
                    end
                end
            end
            if px(i)<0      % min_{px<0} { min(alphaXi) }
                if alphaMinus>alphaXi
                    alphaMinus=alphaXi;
                end
            elseif px(i)>0  % min_{px>0} { min(alphaXi) }
                if alphaPlus>alphaXi
                    alphaPlus=alphaXi;
                end
            end
            %-------------------
            % s(sigma,alphaS)>=psi_k
            %-------------------
            psi_i=min(nu_k,rho*min(s));
            c1=s(i)-psi_i;
            c2=sp(i); alphaSi=pi/2;
            if c2==0 && c1+spp(i) < 0 % Case 1a
                alphaSi=acos((c1+spp(i))/spp(i));
                if alphaSi < alphaS 
                    alphaS=alphaSi;
                    if alphaS < 0
                        warning('check');
                    end
                end
            end
            if spp(i)==0 && c2>c1 % Case 2a
                alphaSi=asin((c1)/c2);
                if alphaSi < alphaS 
                    alphaS=alphaSi;
                    if alphaS < 0
                        warning('check');
                    end
                end
            end
            if c2>0 && spp(i)>0 && c1+spp(i) < sqrt(c2^2+spp(i)^2) % Case 3a
                alphaSi=asin((c1+spp(i))/sqrt(c2^2+spp(i)^2)) ...
                    -asin((spp(i))/sqrt(c2^2+spp(i)^2));
                if alphaSi < alphaS 
                    alphaS=alphaSi;
                    if alphaS < 0
                        warning('check');
                    end
                end
            end
            if c2>0 && spp(i)<0 && c1+spp(i) < sqrt(c2^2+spp(i)^2) % Case 4a
                alphaSi=asin((c1+spp(i))/sqrt(c2^2+spp(i)^2)) ...
                    +asin((-spp(i))/sqrt(c2^2+spp(i)^2));
                if alphaSi < alphaS 
                    alphaS=alphaSi;
                    if alphaS < 0
                        warning('check');
                    end
                end
            end
            if c2<0 && spp(i)<0 && c1+spp(i) < 0 % Case 5a
                alphaSi=pi-asin(-(c1+spp(i))/sqrt(c2^2+spp(i)^2)) ...
                    -asin((-spp(i))/sqrt(c2^2+spp(i)^2));
                if alphaSi < alphaS 
                    alphaS=alphaSi;
                    if alphaS < 0
                        warning('check');
                    end
                end
            end
            if ps(i)<0      % min_{ps<0} { min(alphaSi) }
                if alphaMinus>alphaSi
                    alphaMinus=alphaSi;
                end
            elseif ps(i)>0  % min_{ps>0} { min(alphaSi) }
                if alphaPlus>alphaSi
                    alphaPlus=alphaSi;
                end
            end
        end %end  of for loop 1:n
        %----------------
        % Algo 6.1 Step 3
        %----------------
        if alphaMinus>alphaPlus
            sigma_lb=Sigma;
        elseif alphaMinus<=alphaPlus
            sigma_ub=Sigma;
        end
    % Algo 6.1 Step 4
    end 
    alpha=min(min(alphaMinus,alphaPlus),min(alphaX,alphaS));%in case ps(i)=0 or px(i)=0

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Algo 5.1 Step 2 u>u(sigma,alphaX)
    %%%%%%%%%%%%%%%%%%%%%%%%
    reduceAlpha=1;
    a1=dx'*ds/n;
    a2=(dx'*qs+ds'*qx)/n;
    a3=(dx'*ps+ds'*px)/n;
    while reduceAlpha==1
        reduceAlpha=0;
        ualpha=mu*(1-sin(alpha))+(1-cos(alpha))*((mu-a3*sin(alpha))*Sigma...
            -a1*(1-cos(alpha))-a2*sin(alpha));
        if mu<=ualpha
            alpha=alpha/2;
            ualpha=mu*(1-sin(alpha))+(1-cos(alpha))*((mu-a3*sin(alpha))*Sigma...
                -a1*(1-cos(alpha))-a2*sin(alpha));
            reduceAlpha=1;
        end
    end
    % rescale alpha_k
%     if alpha>maxAlpha
%         alpha=maxAlpha;
%     end
    alpha=min(0.999*alpha,maxAlpha2);
    %%%%%%%%%%%%%%%%%%%%%%
    % Algo 5.1 Step 3 
    %%%%%%%%%%%%%%%%%%%%%%
    ddx=px*Sigma+qx;
    ddl=pl*Sigma+ql;
    dds=ps*Sigma+qs;
    normRbk=norm(A*x-b); normRck=norm(A'*lambda+s-c);
    lambda=lambda-dl*sin(alpha)+ddl*(1-cos(alpha));
    s=s-ds*sin(alpha)+dds*(1-cos(alpha));
    x=x-dx*sin(alpha)+ddx*(1-cos(alpha));
    rB=A*x-b;
    rC=A'*lambda+s-c;
    mu=s'*x/n
    nu_k=nu_k*(1-sin(alpha));       
    %%%%%%%%%%%%%%%%%%%%%%
    % Step 4 
    %%%%%%%%%%%%%%%%%%%%%% 
    kk=kk+1
%     if kk==2
%         warning('check')
%     end
    if alpha <10^(-8) || mu<10^(-8)
        break;
    end
    if (10*normRbk<norm(rB) && norm(rB)>10^(-6)) ...
            || (10*normRck<norm(rC) && norm(rC)>10^(-6))
%         [normRbk norm(rB) normRck norm(rC)]
        x=xOld; s=sOld; mu=s'*x/n;
        break;
    end
    if d==1
        if min(x)<10^(-6) % && mu>10^(-7)
            singleEarly=1;
        end
    end
    sOld=s; xOld=x; lOld=lambda;
end
obj=c'*x+csum;
infe=norm(A*x-b);

% recover back to original x vector
if d==0
    for i=1:length(x_indx)
        x_final(x_indx(i))=x(i);
    end
    for ii=length(x_tobe):-1:1
        idxp=find(A_post(ii,:));
        x_final(x_tobe(ii))=b_post(ii);
        for jj=1:length(idxp)
            if idxp(jj)~=x_tobe(ii)
                x_final(x_tobe(ii))=x_final(x_tobe(ii))- ...
                    A_post(ii,idxp(jj))*x_final(idxp(jj));
            end
        end
    end
    obj=c_orig'*x_final; x=x_final;
else
    obj=c'*x+csum;
end
% finish the recovery



% this function computes the smallest and biggist absolute nonzeros
% which is used to decide if we need to scale the matrix A b c
function [minA,maxA]=calRatioCondition(A)
m=length(A(:,1));
% rowMin=zeros(m,1); rowMax=zeros(m,1);
[i,j,v]=find(A);
if ~isempty(v)
    minA=min(abs(v));
    maxA=max(abs(v));
else
    error('A is a matrix of zero')
end




function [A,b,c,x,s,lambda,rB,rC,lOld,xOld,sOld] = ...
    makeAfull(A,b,c,x,s,lambda,rB,rC,lOld,xOld,sOld)
if exist('lambda','var') == 0
    indlambda = 0;
else
    indlambda = 1;
end
[m,n]=size(A);
indA=sparse(m,n); deA=A; deb=b; indb=zeros(m,1); 
if indlambda==1
    ilambda=zeros(m,1); delambda=lambda; 
    irB=zeros(m,1); derB=rB;
    iold=zeros(m,1); dold=lOld;
end
m1=1; m2=m; sj=1;
% separate independent rows indA and potential dependent row deA
% using column singleton
while sj<n && m1<=m && m2>0
    if nnz(deA(:,sj))==1
        % swap rows of A and make relative changes to b,lambda,rB,
        iw=find(deA(:,sj));
        indA(m1,:)= deA(iw,:); 
        indb(m1)=deb(iw); % align b with A
        if indlambda==1
            ilambda(m1)=delambda(iw); delambda(iw)=[]; % align lambda with A
            irB(m1)=derB(iw); derB(iw)=[]; % align rB with A
            iold(m1)=dold(iw); dold(iw)=[];
        end
        deA(iw,:)=[]; deb(iw)=[]; 
        %swap columns m1 and sj of A and make relative changes to
        %c,s,rC
        if sj~=m1
            indA(:,[m1,sj])=indA(:,[sj,m1]); 
            deA(:,[m1,sj])=deA(:,[sj,m1]);
            c([m1,sj])=c([sj,m1]);
        end
        if indlambda==1
            s([m1,sj])=s([sj,m1]); sOld([m1,sj])=sOld([m1,sj]);
            x([m1,sj])=x([sj,m1]); xOld([m1,sj])=xOld([m1,sj]);
            rC([m1,sj])=rC([sj,m1]);
        end
        sj=m1; m1=m1+1; m2=m2-1; 
    end
    sj=sj+1;
end
if m1-1==m % if A is full rank, b=indb, rB=irB, and lambda=ilambda
    b=indb;
    if indlambda==1
        rB=irB;
        lambda=ilambda; lOld=iold;
    end
else
    %Using Markowitz's pivot rule and pivoting to remove dependent rows in deA
%     deA=deA(1:m2,m1:n); 
%     indA(m1:m,m1:n)=deA; b=[indb; deb];
    indA(m1:m,:)=deA; b=[indb(1:m1-1); deb]; 
    if indlambda==1
        lambda=[ilambda(1:m1-1); delambda];  iOld=[iold(1:m1-1); dold];
        rB=[irB(1:m1-1); derB];
    end
    clear deA deb delambda
    m2=m-m1+1; % number of rows of deA
    while m1<=m 
        minobj=10^10;
        % Markowitz's pivot rule
        for i=m1:m
            ridx=find(abs(indA(i,:))>10^(-8)); 
            % avoid extract small pivot in problem bandm
            rnum=length(ridx);
            if rnum==0  % ith row of A is zero 
                indA(i,:)=[];
                b(i)=[]; 
                if indlambda==1
                    lambda(i)=[]; rB(i)=[]; lOld(i)=[];
                end
                m=m-1;
                break;
            end
            for j=1:rnum % ith row of A is nonzero
    %             cidx=find(indA(m1-1:m,ridx(j)));
    % %%%%%%%%%%%%
                cidx=find(abs(indA(:,ridx(j)))>10^(-8)); 
                % avoid extract small pivot in problem bandm
                rcNum=(rnum-1)*length(cidx);
                if rcNum<minobj
                    minobj=rcNum;
                    minR=i; minRx=ridx;  minCx=cidx; minC=ridx(j);
                end
            end
        end
        if rnum~=0
            % Gauss illimination
            for i=1:length(minCx)
                if minR~=minCx(i)
                    coe=-indA(minCx(i),minC)/indA(minR,minC);
                    b(minCx(i))=b(minCx(i))+coe*b(minR);
                    if indlambda==1
                        rB(minCx(i))=rB(minCx(i))+coe*rB(minR);
                        lambda(minR)=lambda(minR)-coe*lambda(minCx(i));
                        lOld(minR)=lOld(minR)-coe*lOld(minCx(i));
                    end
                    for j=1:length(minRx)
                        indA(minCx(i),minRx(j))=indA(minCx(i),minRx(j)) ...
                            +coe*indA(minR,minRx(j));
                    end
                end
            end
            % column minC is a column singleton
            % Swap rows and columns
            if minR~=m1
                indA([m1,minR],:)=indA([minR,m1],:); %swap rows
                b([m1,minR])=b([minR,m1]);
                if indlambda==1
                    lambda([m1,minR])=lambda([minR,m1]); 
                    lOld([m1,minR])=lOld([minR,m1]);
                    rB([m1,minR])=rB([minR,m1]);
                end
            end
            if m1~=minC
                indA(:,[m1,minC])=indA(:,[minC,m1]); %swap column
                c([m1,minC])=c([minC,m1]);
                if indlambda==1
                    x([m1,minC])=x([minC,m1]); xOld([m1,minC])=xOld([minC,m1]);
                    s([m1,minC])=s([minC,m1]); sOld([m1,minC])=sOld([minC,m1]);
                    rC([m1,minC])=rC([minC,m1]);
                end
            end
            m1=m1+1;
        end
    end
end
% [m,n]=size(indA)
A=indA; 