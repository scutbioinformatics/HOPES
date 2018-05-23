function [X1,X2,X3]=Synthetic_data1(mean,sd);

importdata('..\simulation\rna_sim');
X1=ans.data;
importdata('..\simulation\mirna_sim');
X2=ans.data;
importdata('..\simulation\meth_sim');
X3=ans.data;

%X1 = Standard_Normalization(X1);
%X2 = Standard_Normalization(X2);
%X3 = Standard_Normalization(X3);

XS1=zeros(503,200);
XS1(:,1:50)=2;
XS1(:,151:200)=6;
XS1=XS1+normrnd(0,1,503,200);
X1=reconstruction(X1,XS1);
%X1=Standard_Normalization(X1);
X1=X1+5000*normrnd(mean,sd,503,200);
%noise=normrnd(mean,sd,1000,200);
%X1=[X1;noise];


XS2=zeros(2541,200);

XS2(:,51:100)=4;
XS2(:,151:200)=6;
XS2=XS2+normrnd(0,1,2541,200);
X2=reconstruction(X2,XS2);
%noise=normrnd(mean,sd,5000,200);
%X2=[X2;noise];
%X2=Standard_Normalization(X2);
X2=X2+3.7*normrnd(mean,sd,2541,200);


XS3=zeros(936,200);

XS3(:,1:100)=4;
XS3(:,101:150)=2;
XS3(:,151:200)=6;
XS3=XS3+normrnd(0,1,936,200);
X3=reconstruction(X3,XS3);
%noise=normrnd(mean,sd,2000,200);
%X3=[X3;noise];
X3=X3+0.37*normrnd(mean,sd,936,200);
%X3=Standard_Normalization(X3);


end

function Y = reconstruction(X,XS)
[U,D,V] = svd(X(:,1:200));
%D(2:size(D,2),2:size(D,2))=diag(repmat(trace(D(2:size(D,2),2:size(D,2))/size(D,2)-1),1,size(D,2)-1));
%D(2:size(D,2),2:size(D,2))=0;
[US,DS,VS] = svd(XS);
Y=U*D*VS';
clear U;
clear D;
clear V;
clear US;
clear DS;
clear VS;
end