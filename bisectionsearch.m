function [ Pdcmax ] = bisectionsearch( Pcmax,Pdmax,Pkdmin,alphak,gamma0d,alphamk,noise,p0,epsilon )
%���ַ����Ҹ���λ��
mid=0;
left=Pkdmin;
right=Pdmax;
fleft=(alphak*left/(gamma0d*alphamk)*((exp((-gamma0d*noise)/(left*alphak))/(1-p0))-1))-Pcmax;
fright=(alphak*right/(gamma0d*alphamk)*((exp((-gamma0d*noise)/(right*alphak))/(1-p0))-1))-Pcmax;
fmid=(alphak*((right+left)/2)/(gamma0d*alphamk)*((exp((-gamma0d*noise)/(((right+left)/2)*alphak))/(1-p0))-1))-Pcmax;
if fleft*fright>0
    Pdcmax=Pdmax;
    return;
else%��������������࣬Ӧ�����ҵ�������������Ȼ��0����
while fmid>=epsilon
    mid=(right+left)/2;
    fmid=(alphak*mid/(gamma0d*alphamk)*((exp((-gamma0d*noise)/(mid*alphak))/(1-p0))-1))-Pcmax;
    if fleft*fmid>0
        fleft=fmid;
        left=mid;
    else
        fright=fmid;
        right=mid;
    end
end
end
if mid==0;
    Pdcmax=Pdmax;
else
    Pdcmax=mid;
end
end

