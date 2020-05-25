clc;
close all;
clear all;
load city_31
C=city_31;
N=size(C,1);%31个城市
Dis=zeros(N);%城市间距离矩阵
for i=1:N
    for j=1:N
        Dis(i,j)=((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5;
    end
end
%%参数设置
NP=250;%初始种群规模
GE=10000;%最大迭代次数
f=zeros(NP,N);%初始种群，NP*31
F=[];%中间更新种群存储
%%随机生成初始种群
for i=1:NP
    f(i,:)=randperm(N);
end
lens=zeros(NP,1);%存储路径长度200*1
R=f(1,:);%存储最优路径1*31
fitness=zeros(NP,1);%存储适应度值
%%%%%%%循环遗传算法，最大迭代次数为GE%%%
gen=0;
while gen<GE
    %%计算路径长度
    for i=1:NP
        lens(i,1)=Dis(f(i,N),f(i,1));
        for j=1:N-1
          lens(i,1)=lens(i,1)+Dis(f(i,j),f(i,j+1));
        end
    end
    maxlen=max(lens);
    minlen=min(lens);
    r=find(lens==minlen);
    R=f(r(1,1),:);%更新最优路径
    %%%%%%变化适应度函数%%%%
    for i=1:length(lens)
        fitness(i,1)=(1-((lens(i,1)-minlen)/(maxlen-minlen+0.001)));
    end
    %%%%%%选择算子%%%%%%%%%%%%%%%
    nn=0;
    for i=1:NP
        if fitness(i,1)>=rand
            nn=nn+1;
            F(nn,:)=f(i,:);
        end
    end
    [aa,bb]=size(F);% *31
    while aa<NP
        nnper=randperm(nn);
        A=F(nnper(1),:);
        B=F(nnper(2),:);
    %%%%%%交叉操作%%%%%%%%%%%%%%%
        W=ceil(N/10);%%交叉点个数
        p=unidrnd(N-W+1);
        for i=1:W
            x=find(A==B(p+i-1));%找A中和B待变异范围内元素相等的位置
            y=find(B==A(p+i-1));
            t=A(p+i-1);
            A(p+i-1)=B(p+i-1);
            B(p+i-1)=t;
            t=A(x);
            A(x)=B(y);
            B(y)=t;
        end
    %%%%%%%%%%%变异操作%%%%%%%%%%%
        p1=floor(1+N*rand());
        p2=floor(1+N*rand());
        while p1==p2
             p1=floor(1+N*rand());
             p2=floor(1+N*rand());
        end
        t1=A(p1);
        A(p1)=A(p2);
        A(p2)=t1;
        t1=B(p1);
        B(p1)=B(p2);
        B(p2)=t1;
        F=[F;A;B];
        [aa,bb]=size(F);
    end
    if aa>NP
        F=F(1:NP,:);
    end
    f=F;
    f(1,:)=R;
    clear F;
    gen=gen+1;
    Rlength(gen)=minlen;
end
figure
for i=1:N-1
    plot([C(R(i),1),C(R(i+1),1)],[C(R(i),2),C(R(i+1),2)],'bo-');
    hold on
end
 plot([C(R(N),1),C(R(1),1)],[C(R(N),2),C(R(1),2)],'ro-');
 title(['优化最短距离：',num2str(minlen)]);

 figure
 plot(Rlength)
 xlabel('迭代次数');
 ylabel('目标函数值');
 title('适应度进化值');