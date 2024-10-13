%% �������Ч�����봿�� ����ľ���Ϊ CP���㷨������ʵ�����õ������ݵĽ���
%Ci �㷨����õ���ÿ����������
function [Entropy Purity]=EnAndPur(CP,Ci)
%�õ�����ֵ
[rn cn]=size(CP);
%% ������
%������� precision
for i=1:rn
    for j=1:cn
     precision(i,j)=CP(i,j)/Ci(1,i);    
    end
end
%����ei(i,j)
for i=1:rn
    for j=1:cn
     ei(i,j)=precision(i,j)*log2(precision(i,j));    
    end
end
%
%����ei_sum
for i=1:rn
    ei_sum(i)=-nansum(ei(i,:));
end
%����mi*ei_sum(i)
for j=1:cn
    mmi(j)=Ci(1,j)*ei_sum(j);
end
%����entropy
Entropy=nansum(mmi)/nansum(Ci);
%% ���㴿��Purity
%�ҳ�����һ��
for i=1:rn
     pr_max(i)=max(precision(i,:));    
end
%�����������
for j=1:cn
    nni(j)=Ci(1,j)*pr_max(j);
end
Purity=nansum(nni)/nansum(Ci);
end
