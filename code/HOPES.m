function [W]=HOPES(Wall,K,alpha,beta)
if nargin < 2
    K = 20;
end
if nargin < 3
    alpha = 0.5;
end

if nargin < 4
    beta = 0.25;
end

C = length(Wall);
[m,n]=size(Wall{1});
rho=0.5;

%normalize
for i = 1 : C
    Wall{i} = Wall{i}./repmat(2*(sum(Wall{i},2)-diag(Wall{i})),1,n);
    Wall{i} = (Wall{i} + Wall{i}')/2;
end

%LocalAffinity 
for i = 1 : C
    newW{i} = FindDominateSet(Wall{i},round(K));
end
% initialize Omega and W 
for i = 1 : C
    Omega{i}=zeros(m,n);
end
Wsum=zeros(m,n);
for i = 1 : C
    Wsum = Wsum + Wall{i};
end
P=Wsum/C;

ITER=0;
dual_norm=0;
%last=0;
    for i = 1 : C
        dual_norm = dual_norm+(norm(Wall{i}-P))/C;
    end
ez = 1;
 while(max(dual_norm/rho^2,ez)>0.000001)
    % update W1,...WC
    for i = 1 : C        
        %Wi-SiWj
        %Wall{i}=(2*newW{i}*((Wsum-Wall{i})/(C-1))+rho*P-rho*Omega{i})/(2+rho);
        
        %Wi-SiWjSi^T
        %   Wall{i}=(2*newW{i}*((Wsum-Wall{i})/(C-1))*newW{i}'+rho*P-rho*Omega{i})/(2+rho);
        
        %Wi-SiWjSj^T
        %tmp=zeros(n);
        %for j=1:C
        %   tmp=tmp+(2*newW{i}*((Wsum-Wall{i})/(C-1))*newW{j}');
        %end
        %Wall{i}=(tmp+rho*P-rho*Omega{i})/(2*C+rho);
        
        %W-Si+W-SiW+W-SiWSj^T
        tmp=zeros(n);
        for j=1:C
           tmp=tmp+(2*beta*newW{i}*((Wsum-Wall{i})/(C-1))*newW{j}');
        end
        X=zeros(n);
        X(newW{i}~=0)=1;
        Wall{i}=(2*newW{i}.*X+2*alpha*newW{i}*((Wsum-Wall{i})/(C-1))+tmp+rho*P-rho*Omega{i})./(2.*X.*X+(2*alpha+2*beta*C+rho)*ones(n));
    end
    
%    for i = 1 : C
%        newW{i} = FindDominateSet(Wall{i},round(K));
%    end

    % update W
    Wsum=zeros(m,n);
    for i = 1 : C
        Wsum = Wsum + Wall{i};
     end
    prev_P=P;
    P=Wsum/C;    
    %update Omega
    for i = 1 : C
        Omega{i}=Omega{i}+Wall{i}-P;
    end
    
    temp=0;
    %for i = 1 : C
        %X=zeros(n);
        %X(newW{i}~=0)=1;
        %temp = temp+(norm(P.*X-newW{i},'fro'))^2+alpha*(norm(P-newW{i}*P,'fro'))^2;
       % for j=1:C
        %    temp=temp+beta*(norm(P-newW{i}*P*newW{j}','fro'))^2;
       % end
    %end
    
    %dual_norm(ITER+1) = norm(last-temp);
    %last=temp;
    temp=0;
    for i = 1 : C
        temp = temp+(norm(Wall{i}-P))/C;
    end
    dual_norm(ITER+1) = temp;
    ez(ITER+1) = norm(prev_P-P);
    if dual_norm(ITER+1)>100*ez(ITER+1)
        rho = rho*0.8;
    elseif (ez(ITER+1)>.1*dual_norm(ITER+1) )
        rho = 10*rho;
    end
    ITER=ITER+1;
    %fprintf('iter %d, %f,%f \n',ITER,ez(ITER),dual_norm(ITER));
   % fprintf('iter %d, %f,%f \n',ITER,norm(Wall{1}-P,1),rho);
end

W=P;
W=W./repmat(sum(W,2),1,n);
% W = (W +W'+eye(n))/2;
end

function newW = FindDominateSet(W,K)
[m,n]=size(W);
[YW,IW1] = sort(W,2,'descend');
clear YW;
newW=zeros(m,n);
temp=repmat((1:n)',1,K);
I1=(IW1(:,1:K)-1)*m+temp;
%I1 KNN Index
newW(I1(:))=W(I1(:));
newW=newW./repmat(sum(newW,2),1,n);
%newW=newW.*2;
clear IW1;
clear IW2;
clear temp;
end