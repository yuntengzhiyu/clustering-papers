function ROGC(dataset, num_bases,k,alpha,beta,feareal)

dir = 'datasets/';
dataset = 'glass_uni';
load([dir dataset])
fea=X;
gnd=Y;
if feareal
    fea = NormalizeFea(fea);
end
X=fea';
nclass = length(unique(gnd));
par.num_bases = num_bases;
par.alpha=alpha;
par.beta = beta;
iters=10;
par.num_iters = 100;
warning('off', 'all');
for ave=1:10
try
    num = size(X,2);
    distX = L2_distance_1(X,X);
    [distX1, idx] = sort(distX,2);
    W0 = zeros(num);
    rr = zeros(num,1);
    for i = 1:num
        di = distX1(i,2:k+2);
        rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
        id = idx(i,2:k+2);
        W0(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
    end;
    t=1;
    while t < iters+1
        if t ==1
            [~, S, ~,~,fobj] =GraphSC(X, W0, par.num_bases, par.alpha, par.beta, par.num_iters);
            [la, W1, ~,r] = CAN(S, nclass,k);
            obj(1,t)=fobj+r*(norm(W0,'fro')^2);
            fprintf(('fobj= %f\n'),  obj(1,t));
        else
            [~, S, ~,~,fobj]=GraphSC(X, W1, par.num_bases, par.alpha, par.beta, par.num_iters);
            [la, W1, ~,r] = CAN(S, nclass,k);
            obj(1,t)=fobj+r*(norm(W1,'fro')^2);
            fprintf(('fobj= %f\n'),  obj(1,t));
        end
        t = t + 1;
    end
    result_BestCAN=ClusteringMeasure(gnd,la);
    Final(ave,(1:3))=result_BestCAN;
    disp(['BestCAN. ACC: ',num2str(result_BestCAN(1))]);
    disp(['BestCAN. MIhat: ',num2str(result_BestCAN(2))]);
    disp(['BestCAN. Purity: ',num2str(result_BestCAN(3))]);
    fid = fopen(['Eachresult/',num2str(dataset),'_Each_nBases_',num2str(par.num_bases),'_k_',num2str(k),'_alpha_',num2str(par.alpha),'_beta_',num2str(par.beta),'_fea_',num2str(feareal),'_num_',num2str(ave),'.txt'],'w');
    fprintf(fid,'%d %d %f %f %d %f %f %f \n',par.num_bases,k,par.alpha,par.beta,feareal,result_BestCAN(1),result_BestCAN(2),result_BestCAN(3));
    fclose(fid);
    catch
    result_BestCAN=[0,0,0];
    Final(ave,(1:3))=result_BestCAN;
    disp(['BestCAN. ACC: ',num2str(result_BestCAN(1))]);
    disp(['BestCAN. MIhat: ',num2str(result_BestCAN(2))]);
    disp(['BestCAN. Purity: ',num2str(result_BestCAN(3))]);
    fid = fopen(['Eachresult/',num2str(dataset),'_Each_nBases_',num2str(par.num_bases),'_k_',num2str(k),'_alpha_',num2str(par.alpha),'_beta_',num2str(par.beta),'_fea_',num2str(feareal),'_num_',num2str(ave),'.txt'],'w');
    fprintf(fid,'%d %d %f %f %d %f %f %f \n',par.num_bases,k,par.alpha,par.beta,feareal,result_BestCAN(1),result_BestCAN(2),result_BestCAN(3));
    fclose(fid);
    continue;
end
end
a=mean(Final((1:ave),1));
b=mean(Final((1:ave),2));
c=mean(Final((1:ave),3));
fid = fopen(['result/',num2str(dataset),'_nBases_',num2str(par.num_bases),'_k_',num2str(k),'_alpha_',num2str(par.alpha),'_beta_',num2str(par.beta),'_fea_',num2str(feareal),'.txt'],'w');
fprintf(fid,'%d %d %f %f %d %f %f %f \n',par.num_bases,k,par.alpha,par.beta,feareal ,a,b,c);
fclose(fid);
end


