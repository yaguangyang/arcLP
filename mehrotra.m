%
%                       mehrotra.m Version 1.0
%
%                          Yaguang Yang
%
%
%                          May 2016
%
%
% This function finds the optimal solution for standard linear programming
% problems using Mehrotra's interior-point algorithm. 
%
%
% Synopsis:
%            function [x,obj,kk,infe]=mehrotra(A,b,c,d) 
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
%   m1:   the row dimmension of the reduced A matrix after pro-process.
%   n1:   the column dimmension of the reduced A matrix after pro-process.
%


function [x,obj,kk,infe]=mehrotra(A,b,c,d) 

if exist('d','var') == 0
    d = 0;
end
epsilon=0.00000001; 

c_orig=c;
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
    
%     %rule 2: remove linearly dependent rows
%     i=m;
%     while i>1
%         idxi=find(A(i,:)); Ai=A(i,:);
%         for j=i-1:-1:1
%             idxj=find(A(j,:));
%             if length(idxi)==length(idxj)
%                 if idxi==idxj 
%                     coef=A(j,idxj(1))/Ai(idxi(1));
%                     if norm(coef*Ai-A(j,:))<epsilon^2 && ...
%                                coef*b(i)-b(j)<epsilon
%                         A(j,:)=[]; b(j)=[]; %lambda(j)=[];
%                         m=m-1; i=i-1;
%                         prepro=1;
%                     elseif norm(coef*Ai-A(j,:))<epsilon^2 && ...
%                              coef*b(i)-b(j)>epsilon
%                         error('infeasible problem');
%                     end
%                 end
%             end
%         end
%         i=i-1;
%     end
%     disp('Rule 2: removing duplicate rows')
%     disp([num2str(m),'  ', num2str(n)])

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
    
%     %rule 4: remove duplicate columns
%     i=n;
%     while i>1
%         idxi=find(A(:,i)); Ai=A(:,i); ci=c(i);
%         j=i-1;
%         while j>=1
%             idxj=find(A(:,j));
%             if length(idxi)==length(idxj)
%                 if idxi==idxj
%                     if norm(Ai-A(:,j))<epsilon^2    %A(:,i)-A(:,j)=0
%                         if abs(ci-c(j))<epsilon;           %c(i)-c(j)=0
%                             A(:,j)=[]; c(j)=[]; %s(j)=[]; x(i)=[];
%                             j=j-1; i=i-1; prepro=1; n=n-1;
%                         end
%                     end
%                 end
%             end
%             j=j-1;
%         end
%         i=i-1;
%     end
%     disp('Rule 4: remove duplicate columns')
%     disp([num2str(m),'  ', num2str(n)])

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

%     %rule 6: combined columns that form a free variable, then solve it
%     i=n;
%     while i>1
%         j=i-1; skipJ=0;
%         while j>=1
%             if  norm(A(:,i)+A(:,j))<epsilon^2 & abs(c(i)+c(j))<epsilon^2
%                 %combine two singleton columns into one to become a free
%                 %variable, then solve the free variable, then substitute to
%                 %the objective function
%                 if nnz(A(:,i))==0
%                     if c(i)< -epsilon*100
%                         error('the problem is unbounded')
%                     else
%                         A(:,i)=[]; c(i)=[]; 
%                         A(:,j)=[]; c(j)=[];
%                         i=i-2; n=n-2; prepro=1;
%                         skipJ=1;
%                     end
%                 else
%                     A(:,i)=[]; c(i)=[]; %s(i)=[]; x(i)=[];
%                     i=i-1; n=n-1;
%                     firstRow=1;
%                     for k=1:m
%                         if A(k,j)~=0 %find the first row s.t. A(k,i)~=0
%                             if firstRow~=1
%                                 if A(k,j)~=0 & b1~=0
%                                     b(k)=b(k)-A(k,j)*b1/A(k1,j);
%                                 end
%                                 for ii=1:n
%                                     if ii~=j
%                                         if (A(k,j)~=0 & A(k1,ii)~=0)
%                                             A(k,ii)=A(k,ii)- ...
%                                                  A(k,j)*A(k1,ii)/A(k1,j);
%                                         end
%                                     end
%                                 end
%                             else
%                                 b1=b(k); k1=k;
%                                 if c(j)~=0 & b(k) ~=0
%                                     csum=csum+c(j)*b(k)/A(k,j);
%                                 end
%                                 for ii=1:n
%                                     if ii~=j
%                                         if A(k,ii)~=0 & c(j)~=0
%                                             c(ii)=c(ii)-c(j)*A(k,ii)/A(k,j);
%                                         end
%                                     end
%                                 end
%                                 firstRow=0;
%                             end
%                         end
%                     end
%                     A(:,j)=[]; c(j)=[]; %s(j)=[]; x(j)=[];
%                     j=j-1; n=n-1;
%                     A(k1,:)=[]; b(k1)=[]; m=m-1;
%                      prepro=1; skipJ=1;
%                 end
%             end
%             if skipJ==1
%                 break;
%             else
%                 j=j-1;
%             end
%         end
%         i=i-1;
%     end
%     disp('Rule 6: remove free variabls')
%     disp([num2str(m),'  ', num2str(n)])

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
 
%     %rule 8: for b(i)=b(j), and A(i,:)-A(j,:)>=0 (A(i,:)-A(j,:)<=0) remove  
%     %variables xk corresponding to A(j,k)>0 (A(j,k)<0)
%     i=m;
%     while i>1
%         j=i-1;
%         while j>=1
%             Aimj=A(i,:)-A(j,:); 
%             if b(i)-b(j)==0 & min(Aimj)>=0
%                 nnzAi=nnz(A(i,:)); nnzAj=nnz(A(j,:));
%                 idij=find(Aimj>0);
%                 for k=length(idij):-1:1
%                     A(:,idij(k))=[]; c(idij(k))=[]; n=n-1; prepro=1;
%                 end
%                 if nnzAi>nnzAj
%                     A(i,:)=[]; b(i)=[];
%                     m=m-1; %i=i-1; 
%                 else
%                     A(j,:)=[]; b(j)=[];
%                     m=m-1; %j=j-1; i=i-1; 
%                 end
%             elseif (b(i)-b(j)==0 & max(Aimj)<=0)
%                 nnzAi=nnz(A(i,:)); nnzAj=nnz(A(j,:));
%                 idij=find(Aimj<0);
%                 for k=length(idij):-1:1
%                     A(:,idij(k))=[]; c(idij(k))=[]; n=n-1; prepro=1;
%                 end
%                 if nnzAi>nnzAj
%                     A(i,:)=[]; b(i)=[];
%                     m=m-1; %i=i-1; 
%                 else
%                     A(j,:)=[]; b(j)=[];
%                     m=m-1; %j=j-1; i=i-1; 
%                 end  
%             end
%             j=j-1;
%         end
%         i=i-1;
%     end
%     i=m;
%     while i>1
%         j=i-1;
%         while j>=1
%             Aipj=A(i,:)+A(j,:);
%             if b(i)+b(j)==0 & min(Aipj)>=0
%                 nnzAi=nnz(A(i,:)); nnzAj=nnz(A(j,:));
%                 idij=find(Aipj>0);
%                 for k=length(idij):-1:1
%                     A(:,idij(k))=[]; c(idij(k))=[]; n=n-1; prepro=1;
%                 end
%                 if nnzAi>nnzAj
%                     A(i,:)=[]; b(i)=[];
%                     m=m-1; %i=i-1; 
%                 else
%                     A(j,:)=[]; b(j)=[];
%                     m=m-1; %j=j-1; i=i-1; 
%                 end
%             elseif (b(i)+b(j)==0 & max(Aipj)<=0)
%                 nnzAi=nnz(A(i,:)); nnzAj=nnz(A(j,:));
%                 idij=find(Aipj<0);
%                 for k=length(idij):-1:1
%                     A(:,idij(k))=[]; c(idij(k))=[]; n=n-1; prepro=1;
%                 end
%                 if nnzAi>nnzAj
%                     A(i,:)=[]; b(i)=[];
%                     m=m-1; %i=i-1; 
%                 else
%                     A(j,:)=[]; b(j)=[];
%                     m=m-1; %j=j-1; i=i-1; 
%                 end  
%             end
%             j=j-1;
%         end
%         i=i-1;
%     end
%     disp('Rule 8: remove fixed variables defined by two rows')
%     disp([num2str(m),'  ', num2str(n)])
    
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

%     %rule 10: remove fixed variables
%     % two lines are almost the same except one extra element
%     i=m;
%     while i>=1
%         j=i-1;
%         while j>=1
%             tmpRow=A(i,:)-A(j,:); % i and j are both rows
%             if nnz(tmpRow)==1
%                 k=find(tmpRow);
%                 xk=(b(i)-b(j))/tmpRow(k);
%                 if xk>=0
%                     for ell=1:length(b)
%                         if ell~=i & ell~=j
%                             b(ell)=b(ell)-A(ell,k)*xk;
%                         end
%                     end
%                     A(j,:)=[]; b(j)=[]; i=i-1; %lambda(j)=[];
%                     if c(k)~=0 & xk~=0
%                         csum=csum+c(k)*xk;
%                     end
%                     A(:,k)=[]; c(k)=[]; %s(k)=[]; x(k)=[]; 
%                     prepro=1;
%                     m=m-1; n=n-1;
%                 else
%                     error('the problem is infeasible')
%                 end
%             end
%             j=j-1;
%         end
%         i=i-1;
%     end  
%     disp('Rule 10: remove singletons defined by two rows')
%     disp([num2str(m),'  ', num2str(n)])
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

%----------------------------------------
%scaling of A, b, c if A is poorly scaled
%----------------------------------------
%calculate the ratio of max(|aij)|/min(|aij|) for aij~=0
% [minApre,maxApre]=calRatioCondition(A);
% Aratio=maxApre/minApre;
% if Aratio > 10^6
%     Abackup=A; bbackup=b; cbackup=c;
%     % form normal equation (see Curtis and Reid paper)
%     E=sparse(m,n);
%     M=sparse(m,m); p1=zeros(m,1);
%     for i=1:m
%         M(i,i)=nnz(A(i,:));
%         idxj=find(A(i,:));
%         for j=1:length(idxj)
%             p1(i)=p1(i)+log10(abs(A(i,idxj(j))));
%             E(i,idxj(j))=1;
%         end
%     end
%     N=sparse(n,n); p2=zeros(n,1);
%     for i=1:n
%         N(i,i)=nnz(A(:,i));
%         idxj=find(A(:,i));
%         for j=1:length(idxj)
%             p2(i)=p2(i)+log10(abs(A(idxj(j),i)));
%         end
%     end
%     BA=[M E; E' N];
%     pb=[p1; p2];
%     %inverse constantly be used
%     invM=inv(M); invN=inv(N); pConv=trace(M);
%         A=Abackup; b=bbackup; c=cbackup;
%         y0=[invM*p1; zeros(n,1)];
%         if Aratio < 10^7
%             stp=6;
%         elseif Aratio < 10^8
%             stp=9;
%         elseif Aratio < 10^9
%             stp=12;
%         else
%             stp=15;
%         end
%         [f y]=fscaling(A,M,E,N,p1,p2,y0,m,n,stp);
%         sca1=y(1:m);
%         sca2=y(m+1:m+n);
%         sR1=(1./(10.^sca1)); 
%         sC1=(1./(10.^sca2));
%         [A,b,c]=takeScaling(A,b,c,sR1,sC1);
%         [minA,maxA]=calRatioCondition(A)
%         maxA/minA
%         disp('deepest decent rescaling is used')
%         clear pb ss rowMin rowMax I J S sR1 sC1 idxi idxj sca ...
%             sca1 sca2 BA E M N p1 p2 invM invN Abackup bbackup cbackup;
% else 
%     disp('rescaling is not used')
% end

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
end
xini=x; sini=s;


mu=s'*x/n; rB=A*x-b; rC=A'*lambda+s-c;
kk=0; %iteration count
singleEarly=0;
sOld=s; xOld=x; lOld=lambda;
while (norm(rB)/max(1, norm(b)))+(norm(rC)/max(1, norm(c))) ...
        +(abs(c'*x-b'*lambda)/max(norm(b),max(1,norm(c))))>epsilon
    if d==1
        if fullRank==0 || singleEarly==1
            disp('real initial x-i-b s-i-c')
%             [max(x) max(xini) max(b) max(s) max(sini) max(c)]
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
    % Algorithm 3.3 Step 1 
    %%%%%%%%%%%%%%%%%%%%%%%%

    %------------------------------------------------------------------
    % dx dxx dlambda ddlambda ds dss represented by polynomial of sigma
    %------------------------------------------------------------------
    invS=sparse(diag(1./s)); 
    XinvS=sparse(diag(x.*(1./s))); 
    AXinvS=sparse(A*XinvS);
    tmpA=sparse(AXinvS*A');
    tmpb=AXinvS*rC-b;
    dl=tmpA\tmpb;
    ds=rC-A'*dl;
    dx=x-XinvS*ds;

    Alpha_x=1;
    for i=1:n
        if dx(i)>0
            midAlpha=x(i)/dx(i);
            if midAlpha<Alpha_x
                Alpha_x=midAlpha;
            end
        end
    end
    Alpha_s=1;
    for i=1:n
        if ds(i)>0
            midAlpha=s(i)/ds(i);
            if midAlpha<Alpha_s
                Alpha_s=midAlpha;
            end
        end
    end
    mu_aff=((x-dx*Alpha_x)'*(s-ds*Alpha_s))/n;
    
    
    xigma=min((mu_aff/mu)^3,0.5);

    tmpb=A*invS*(2*dx.*ds-xigma*mu*e);
    ddl=tmpA\tmpb;
    dds=-A'*ddl;
    ddx=invS*(xigma*mu*e-x.*dds-2*(dx.*ds)); 

    clear tmpA tmpb Abackup bbackup
 
    maxAlpha=1; minAlpha=0; Alpha_x=minAlpha;
    while (maxAlpha-minAlpha) > 10^(-6)
        midAlpha=(maxAlpha+minAlpha)/2;
        xAlpha=x-dx*midAlpha+ddx*midAlpha;
        if min(xAlpha) <=0 
            maxAlpha=midAlpha;
        else
            Alpha_x=midAlpha;
            minAlpha=midAlpha;
            xNor=xAlpha;
        end
    end 

    c1=1-exp(-(kk+2)); c2=1+exp(-(kk+2));
    Alpha_1=max(c1*Alpha_x,epsilon); Alpha_2=min(c2*Alpha_x,1);
    xAlpha1=x-dx*Alpha_1+ddx*Alpha_1;
    xAlpha2=x-dx*Alpha_2+ddx*Alpha_2;
    if min(xAlpha1)>min(xNor)
        Alpha_x=Alpha_1; %xAlpha=xAlpha1;
    end
    if min(xAlpha2)>min(min(xNor),min(xAlpha1))
        Alpha_x=Alpha_2; %xAlpha=xAlpha2;
    end
    maxAlpha=1; minAlpha=0; Alpha_s=minAlpha;
    while (maxAlpha-minAlpha) > 10^(-6)
        midAlpha=(maxAlpha+minAlpha)/2;
        sAlpha=s-ds*midAlpha+dds*midAlpha;
        if min(sAlpha) <= 0 
            maxAlpha=midAlpha;
        else
            Alpha_s=midAlpha;
            minAlpha=midAlpha;
            sNor=sAlpha;
        end
    end
    
    Alpha_1=max(c1*Alpha_s,epsilon); Alpha_2=min(c2*Alpha_s,1);
    sAlpha1=s-ds*Alpha_1+dds*Alpha_1;
    sAlpha2=s-ds*Alpha_2+dds*Alpha_2;
    if min(sAlpha1)>min(sNor)
        Alpha_s=Alpha_1; %sAlpha=sAlpha1;
    end
    if min(sAlpha2)>min(min(sNor),min(sAlpha1))
        Alpha_s=Alpha_2; %sAlpha=sAlpha2;
    end
    x=x-dx*Alpha_x+ddx*Alpha_x;
    s=s-ds*Alpha_s+dds*Alpha_s;
    lambda=lambda-dl*Alpha_s+ddl*Alpha_s;
    normRbk=norm(A*x-b); normRck=norm(A'*lambda+s-c);
    
    if min(x)<0 || min(s)<0
        error('DEBUG: something is wrong')
    end
        
%     Sp=sparse(diag(s));
%     Xp=sparse(diag(x));

    rB=A*x-b;
    rC=A'*lambda+s-c;
    mu=s'*x/n;
    
    kk=kk+1;
    if Alpha_x <10^(-8) && Alpha_s<10^(-8)
        break;
    end
    if (10*normRbk<norm(rB) && norm(rB)>10^(-6)) || (10*normRck<norm(rC) ...
            && norm(rC)>10^(-6))
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


% this function computes scaling matrices using deepest descent method
% and golden section search
function [f,y]=fscaling(A,M,E,N,p1,p2,y0,m,n,stp)
fy0=fyfunction(A,y0);
fpy0=[M E;E' N]*y0-[p1;p2]; ik=0;
while norm(fpy0)>0.01 && ik<stp
    y=y0-fpy0;
    fy=fyfunction(A,y);
    [fy,y] = myGolden(A,y0,fpy0);
    y0=y;
    fy0=fy;
    fpy0=[M E;E' N]*y0-[p1;p2]; ik=ik+1;
end
f=fy;

% golden section is used for gradient method to find scaling matrices
function [fy,y] = myGolden(A,y0,fpy0)
a=0; b=1;                       % start of interval
iter= 15;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0; 
x1=a+(1-tau)*(b-a);             % computing x values
x2=a+tau*(b-a);
y1=y0-x1*fpy0;
f_x1=fyfunction(A,y1);
y2=y0-x2*fpy0;
f_x2=fyfunction(A,y2);
while ((abs(b-a)>0.001) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        b=x2;
        x2=x1;
        x1=a+(1-tau)*(b-a);
        y1=y0-x1*fpy0;
        f_x1=fyfunction(A,y1);
        y2=y0-x2*fpy0;
        f_x2=fyfunction(A,y2);
    else
        a=x1;
        x1=x2;
        x2=a+tau*(b-a);
        y1=y0-x1*fpy0;
        f_x1=fyfunction(A,y1);
        y2=y0-x2*fpy0;
        f_x2=fyfunction(A,y2);
    end
    k=k+1;
end
if(f_x1<f_x2)
    fy=f_x1;
    y=y1;
else
    fy=f_x2;
    y=y2;
end

% given y, this function computes f(y) used in fscaling function
function fy=fyfunction(A,y)
m=length(A(:,1)); fy=0;
for i=1:m;
    idxj=find(A(i,:));
    for j=1:length(idxj)
        fy=fy+(log10(abs(A(i,idxj(j))))-y(i)-y(m+idxj(j)))^2;
    end
end

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

% given row scaling sR1 and column scaling matrices
% this function scale the matrices A b c
function [A,b,c]=takeScaling(A,b,c,sR1,sC1)
[m,n]=size(A);
for i=1:m
    b(i)=b(i)*sR1(i);
    idxj=find(A(i,:));
    for j=1:length(idxj)
        A(i,idxj(j))=sR1(i)*A(i,idxj(j));
    end
end
for i=1:n
    c(i)=c(i)*sC1(i);
    idxj=find(A(:,i));
    for j=1:length(idxj)
        A(idxj(j),i)=sC1(i)*A(idxj(j),i);
    end
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
