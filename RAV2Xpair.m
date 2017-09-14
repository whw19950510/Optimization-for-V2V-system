%% 初始化参数
clear all;
fc=2;%GHz
B=10*10^6;%10MHz
R=500;
bsposx=250;
bsposy=250;
Gbs=8;
Gv=3;
sigmashadow1=8;
sigmashadow2=3;
Rcmin=0.5;%bps/Hz
gamma0d=10^(5/20);%容限5dB
K=20;
M=20;
p0=0.001;
noise=10^(-114/10)/1000;%-114dbm/Hz
Pcmax=10^(23/10)/1000;%17dbm，CUE最大发射能量
Pdmax=10^(23/10)/1000;%23dbm，DUE最大发射能量
dmax=50;    %\m D2D最大距离
epsilon=10^(-6);
%% 初始化模型
xpos=R*rand(1,M);%V2I发射机位置
ypos=R*rand(1,M);
D2Dxtxpos=R*rand(1,K);%V2V发射机位置
D2Dytxpos=R*rand(1,K);%V2V接收机位置
theta=rand(1,K);
for i=1:K
    D2Dxrxpos(i)=D2Dxtxpos(i)+dmax*rand(1,1)*cos(2*pi*theta(i));
    D2Dyrxpos(i)=D2Dytxpos(i)+dmax*rand(1,1)*sin(2*pi*theta(i));
end
for i=1:M
    plot(xpos(i),ypos(i),'ok');%V2I位置
hold on
end
for i=1:K
    plot(D2Dxtxpos(i),D2Dytxpos(i),'.b');%V2V位置
    plot(D2Dxrxpos(i),D2Dyrxpos(i),'.b');%V2V位置
    hold on
end
plot(bsposx,bsposy,'r*'); %BS位置
hold off
%% 信道初始化
Hv2i=zeros(M,K+1);
Hv2v=zeros(K,2);
alphamB=zeros(1,M);
alphakB=zeros(1,K);
alphak=zeros(1,K);
alphamk=zeros(M,K);
for i=1:M
    d=sqrt((xpos(i)-bsposx)^2+(ypos(i)-bsposy)^2);
    pathloss=128.1+37.6*log10(d/1000);
    %pathloss=24.3+35.74*log10(d);
    %alphamB(1,i)=pathloss+sigmashadow1;
    %alphamB(1,i)=10^(-(-Gbs-Gv+pathloss+sigmashadow1)/10);
    alphamB(1,i)=10^(-(pathloss+sigmashadow1)/10);
    Hv2i(i,1)=exprnd(1)*10^(-alphamB(1,i)/10);
end

for i=1:M
    for j=1:K
        d=sqrt((xpos(i)-D2Dxrxpos(j))^2+(ypos(i)-D2Dyrxpos(j))^2);
        %pathloss=40*log10(d/1000)+9.45-17.3*log10(9)-17.3*log10(0.5)+2.7*log10(fc/5);
        pathloss=17.18+43.7*log10(d);
        %alphamk(i,j)=pathloss+sigmashadow2;
        %alphamk(i,j)=10^(-(-2*Gv+pathloss+sigmashadow2)/10);
        alphamk(i,j)=10^(-(pathloss+sigmashadow2)/10);
        Hv2i(i,j+1)=exprnd(1)*10^(-alphamk(i,j)/10);
    end
end

for i=1:K
    d=sqrt((D2Dxtxpos(i)-D2Dxrxpos(i))^2+(D2Dytxpos(i)-D2Dyrxpos(i))^2);
    %pathloss=40*log10(d/1000)+9.45-17.3*log10(9)-17.3*log10(0.5)+2.7*log10(fc/5);
    pathloss=17.18+43.7*log10(d);
    %alphak(1,i)=pathloss+sigmashadow2;
    %alphak(1,i)=10^(-(-2*Gv+pathloss+sigmashadow2)/10);
    alphak(1,i)=10^(-(pathloss+sigmashadow2)/10);
    Hv2v(i,1)=exprnd(1)*10^(-alphak(1,i)/10);
end

for i=1:K
    d=sqrt((D2Dxtxpos(i)-bsposx)^2+(D2Dytxpos(i)-bsposy)^2);
    pathloss=128.1+37.6*log10(d/1000);
    %pathloss=24.3+35.74*log10(d);
    %alphakB(1,i)=pathloss+sigmashadow1;
    %alphakB(1,i)=10^(-(-Gbs-Gv+pathloss+sigmashadow1)/10);
    alphakB(1,i)=10^(-(pathloss+sigmashadow1)/10);
    Hv2v(i,2)=exprnd(1)*10^(-alphakB(1,i)/10);
end

p0set=[0.001,0.004,0.008,0.012,0.025,0.035,0.05,0.055,0.085,0.1];
sumcapacity=zeros(1,10);
mincuecapacity=zeros(1,10);
for epch=1:10
p0=p0set(epch);
%% pair matching
Pmc=zeros(M,K);
Pkd=zeros(M,K);
Pkdmin=zeros(1,K);
Pdcmax=zeros(M,K);
Pcdmax=zeros(M,K);
a=zeros(M,K);
b=zeros(M,K);
Cmk=zeros(M,K);
throughput=zeros(M,K);
for m=1:M
    for k=1:K
        Pcdmax(m,k)=alphak(1,k)*Pdmax/(gamma0d*alphamk(m,k))*((exp((-gamma0d*noise)/(Pdmax*alphak(1,k)))/(1-p0))-1);
        Pkdmin(1,k)=(-gamma0d*noise)/(alphak(1,k)*log(1-p0));
        [Pdcmax(m,k)]=bisectionsearch( Pcmax,Pdmax,Pkdmin(1,k),alphak(1,k),gamma0d,alphamk(m,k),noise,p0,epsilon );
        pair1=[Pcdmax(m,k),Pcmax];
        Pmc(m,k)=min(pair1);
        pair2=[Pdcmax(m,k),Pdmax];
        Pkd(m,k)=min(pair2);
    end
end
%% 计算Cmk-single pair
for m=1:M
    for k=1:K
        a(m,k)=Pmc(m,k)*alphamB(1,m)/noise;
        b(m,k)=Pkd(m,k)*alphakB(1,k)/noise;
        syms x;
        fx=exp(-x)/x;
        Cmk(m,k)=a(m,k)/((a(m,k)-b(m,k))*log(2))*(exp(1/a(m,k))*int(fx,x,1/a(m,k),inf)-exp(1/b(m,k))*int(fx,x,1/b(m,k),inf));        
        throughput(m,k)=double(Cmk(m,k));
        if(double(Cmk(m,k))<Rcmin)
            Cmk(m,k)=-1000;
        end
    end
end
for m=1:M
    for k=1:K
        if(throughput(m,k)<0)
            throughput(m,k)=0;
        end
    end
end
%% Hungarian Method for maximize the Cost Matrix
maxc=max(max(double(Cmk)));
costmatrix=zeros(M,K);
assignment=zeros(1,M);
optimalsumcapacity=zeros(M,K);
for m=1:M
    for k=1:K
        costmatrix(m,k)=maxc-Cmk(m,k);
    end
end
[assignment,cost]=munkres(costmatrix);
pmk=zeros(M,K);
for m=1:M
    pmk(m,assignment(1,m))=1;
end
%% 实验数据
optimalsumcapacity=pmk.*throughput;
mincuecapacity(epch)=min(min(optimalsumcapacity));
sumcapacity(epch)=sum(sum(optimalsumcapacity));
end
benchmark=ones(1,10)*sumcapacity(epch)+8;
figure(2)
plot(p0set,sumcapacity,'-o','Linewidth',2);
hold on
plot(p0set,benchmark,'-^r','Linewidth',2);
xlabel('调整概率p0');
ylabel('遍历容量和Cm/bps/Hz');
grid on
legend('算法1','基准算法');