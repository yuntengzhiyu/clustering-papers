% 
clear;
classnum = 7;
each_class_num = 30;
viewnum = 5;
load('X.mat');
% L2 normalization to original data
if 0
    for v = 1:viewnum
        X(v) = {X{v}./norm(X{v},2)};
    end
end

% gaussian normalization
if 0
   for v = 1 :viewnum
       for  j = 1 : classnum*each_class_num;
            data{v}(:,j) = ( data{v}(:,j) - mean( data{v}(:,j) ) ) / std( data{v}(:,j) ) ;     
       end 
   end
end

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
%  vv = rand(1,viewnum);
%  alpha = vv/sum(vv);
NITER = 100;
lambda = 1;

for iter = 1:NITER
    
    % fix alpha, update S
    if iter ==1
       [y, S] = CLR(alpha,A,classnum,lambda);
    else
       [y, S] = CLR(alpha,A,classnum,lambda,S0);
    end
    result_(iter,:) = ClusteringMeasure(groundtruth,y); %ACC NMI Purity
     
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
myresult = ClusteringMeasure(groundtruth,y); %ACC NMI Purity

