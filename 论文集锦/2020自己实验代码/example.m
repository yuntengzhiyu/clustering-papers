clear;	
load('dermatology_uni.mat');
nClass = length(unique(Y));
X = NormalizeFea(X);

W = constructW(X);

num = 0;
c = 6;

while(num < 30)  %ֱ���������
    nBasis = 20;
    alpha = 0.01;
    beta = 0.1;
    nIters = 15;
    %���ݼ������W�����Ը���S,B
    [B, S, stat] = GraphSC(X', W, nBasis, 2 * alpha, beta, nIters); 
    %����F��W
    [la, A, evs2] = CAN(S, c, W, alpha);
    W = A;
    num = num + 1;
end


result = ClusteringMeasure(Y,la);
fprintf('%.4f %.4f %.4f\n',result(1),result(2),result(3));


