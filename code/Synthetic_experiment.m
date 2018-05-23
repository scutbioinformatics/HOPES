function accuracy_our=Synthetic_experiment(mean,sd,alpha,beta);
if nargin < 2
    sd = 2;
end
if nargin < 3
    alpha = 0.5;
end

if nargin < 4
    beta = 0.25;
end
K = 15;%neighbor number

%ground_truth
numOfModules = 4;
n = 200;
ground_truth = [];
for k = 1 : numOfModules
    ground_truth = [ground_truth; k * ones(n / numOfModules, 1)];
end
%data construction
[X1,X2,X3]=Synthetic_data1(mean,sd);

Data1 = Standard_Normalization(X1');
Data2 = Standard_Normalization(X2');
Data3 = Standard_Normalization(X3');

%Data1 =X1';
%Data2 =X2';
%Data3 =X3';

Dist1 = dist2(Data1,Data1);
Dist2 = dist2(Data2,Data2);
Dist3 = dist2(Data3,Data3);

W1 = affinityMatrix(Dist1, K);
W2 = affinityMatrix(Dist2, K);
W3 = affinityMatrix(Dist3, K);

group1 = SpectralClustering(W1,3);
displayClusters(W1,group1,0);

group2 = SpectralClustering(W2,3);
displayClusters(W2,group2,0);

group3 = SpectralClustering(W3,3);
displayClusters(W3,group3,0);

W_our = HOPES({W1,W2,W3},K,alpha,beta);
group_our = SpectralClustering(W_our,4);
    IDX=zeros(10,size(group_our,1));
    for i=1:10
        try
            IDX(i,:)=SpectralClustering(W_our,4);
        catch
            IDX(i,:)=group_our;
        end
    end
    U = {'U_H','std',[]};   
    K = 4; % number of clusters for consensus clustering
    r = 10;
    w = ones(r,1); % the weight of each partitioning
    rep = 10; % the number of ECC runs
    maxIter = 40;
    minThres = 1e-5;
    utilFlag = 0;
    [pi_sumbest,pi_index,pi_converge,pi_utility,t] = RunECC(IDX',K,U,w,rep,maxIter,minThres,utilFlag); % run KCC for consensus clustering

accuracy_our = nmi(pi_index, ground_truth);
displayClusters(W_our,group_our,0);








