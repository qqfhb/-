clc;
close all;
clear all;
load city_31
C=city_31;
N=size(C,1);%31������
Dis=zeros(N);%���м�������
for i=1:N
    for j=1:N
        Dis(i,j)=((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5;
    end
end
%%��������
NP=250;%��ʼ��Ⱥ��ģ
GE=10000;%����������
f=zeros(NP,N);%��ʼ��Ⱥ��NP*31
F=[];%�м������Ⱥ�洢
%%������ɳ�ʼ��Ⱥ
for i=1:NP
    f(i,:)=randperm(N);
end
lens=zeros(NP,1);%�洢·������200*1
R=f(1,:);%�洢����·��1*31
fitness=zeros(NP,1);%�洢��Ӧ��ֵ
%%%%%%%ѭ���Ŵ��㷨������������ΪGE%%%
gen=0;
while gen<GE
    %%����·������
    for i=1:NP
        lens(i,1)=Dis(f(i,N),f(i,1));
        for j=1:N-1
          lens(i,1)=lens(i,1)+Dis(f(i,j),f(i,j+1));
        end
    end
    maxlen=max(lens);
    minlen=min(lens);
    r=find(lens==minlen);
    R=f(r(1,1),:);%��������·��
    %%%%%%�仯��Ӧ�Ⱥ���%%%%
    for i=1:length(lens)
        fitness(i,1)=(1-((lens(i,1)-minlen)/(maxlen-minlen+0.001)));
    end
    %%%%%%ѡ������%%%%%%%%%%%%%%%
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
    %%%%%%�������%%%%%%%%%%%%%%%
        W=ceil(N/10);%%��������
        p=unidrnd(N-W+1);
        for i=1:W
            x=find(A==B(p+i-1));%��A�к�B�����췶Χ��Ԫ����ȵ�λ��
            y=find(B==A(p+i-1));
            t=A(p+i-1);
            A(p+i-1)=B(p+i-1);
            B(p+i-1)=t;
            t=A(x);
            A(x)=B(y);
            B(y)=t;
        end
    %%%%%%%%%%%�������%%%%%%%%%%%
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
 title(['�Ż���̾��룺',num2str(minlen)]);

 figure
 plot(Rlength)
 xlabel('��������');
 ylabel('Ŀ�꺯��ֵ');
 title('��Ӧ�Ƚ���ֵ');