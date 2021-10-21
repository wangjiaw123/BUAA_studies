clear,clc  
A=100*rand(3,100); %����100���ռ��е������
[m1,n1]=size(A);
for i=1:n1
    for j=1:n1
        cc=0;
        for k=1:m1
            cc=cc+(A(k,i)-A(k,j))^2;    %�������������ƽ����
        end
        a(i,j)=sqrt(cc);    %������100����֮��ľ���
    end
end


for i=1:n1
    for j=1:n1
        D(i,j)=mydijkstra(a,i,j);  %����dijkstra�㷨�������������֮��ľ������D
        %D(i,j)=myfloyd(a,i,j); %����floyd�㷨�������������֮��ľ������D
    end
end


% r=25;%���õ�ά�ռ��ά��
m=length(D);
aa=0;
B=zeros(m);
for i=1:m
    for j=1:m
        aa=aa+D(i,j)^2;%������ڲ�Ԫ�ص�ƽ����
    end
end
dist_tdd=(1/m^2)*aa;
for i=1:m
    for j=1:m
        dist_td=(1/m)*sum(D(i,:).^2);%���к�
        dist_dj=(1/m)*sum(D(:,j).^2);%���к�
        B(i,j)=-0.5*(D(i,j)^2-dist_td-dist_dj+dist_tdd);%����MDS�㷨����B
    end
end
[V,Lamda]=eig(B)
[LamdaPX,I]=sort(diag(Lamda),'descend')

r=1;
pa=LamdaPX(r);
while sum(pa)/sum(LamdaPX)<1-10^(-15)     %��������ֵ��ռ�ı�ѡ���r
      r=r+1;
      pa=[pa LamdaPX(r)];
end
Vs=V(:,I(1:r));
Lamdas=Lamda(I(1:r),I(1:r));
Z=(Lamdas.^1/2)*Vs';