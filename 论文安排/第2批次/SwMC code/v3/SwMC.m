clear;
classnum = 7;
each_class_num = 30;
viewnum = 5;
load('X.mat');

for v = 1:viewnum
    A(v) = {constructW_PKN(X{1,v},10)};
end

% groundtruth
groundtruth = zeros(size(A{1},1),1); % column vector
for c = 1:classnum
    for cnt = 1:each_class_num
        groundtruth((c-1)*each_class_num+cnt,1) = c;
    end
end
%% outer loop
alpha = 1/viewnum*ones(1,viewnum);
NITER = 100;
lambda = 0;

for iter = 1:NITER
    
    % fix alpha, update S
    if iter ==1
       [y, S] = CLR(alpha,A,classnum,lambda);
    else
       [y, S] = CLR(alpha,A,classnum,lambda,S0);
    end
    
    % fix S, uodate alpha
    for v = 1:viewnum
        alpha(1,v) = 0.5/norm(S-A{v},'fro');
    end
    S0 = S;
    % calculate obj
    obj = 0;
    for v = 1:viewnum
        obj = obj+norm(S-A{v},'fro');
    end
    Obj(iter) = obj;
    if (iter>1 && abs(Obj(iter-1)-Obj(iter)) < 10^-10)
        break;
    end  
end
    lambda = 10;
    [y, S] = CLR(alpha,A,classnum,lambda,S0);
    myresult = ClusteringMeasure(groundtruth,y); %ACC NMI Purity

