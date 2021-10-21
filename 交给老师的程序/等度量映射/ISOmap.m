clear,clc  
A=100*rand(3,100); %生成100个空间中点的坐标
[m1,n1]=size(A);
for i=1:n1
    for j=1:n1
        cc=0;
        for k=1:m1
            cc=cc+(A(k,i)-A(k,j))^2;    %求任意两点距离平方和
        end
        a(i,j)=sqrt(cc);    %计算这100个点之间的距离
    end
end


for i=1:n1
    for j=1:n1
        D(i,j)=mydijkstra(a,i,j);  %调用dijkstra算法计算出任意两点之间的距离矩阵D
        %D(i,j)=myfloyd(a,i,j); %调用floyd算法计算出任意两点之间的距离矩阵D
    end
end


% r=25;%设置低维空间的维度
m=length(D);
aa=0;
B=zeros(m);
for i=1:m
    for j=1:m
        aa=aa+D(i,j)^2;%求矩阵内部元素的平方和
    end
end
dist_tdd=(1/m^2)*aa;
for i=1:m
    for j=1:m
        dist_td=(1/m)*sum(D(i,:).^2);%求行和
        dist_dj=(1/m)*sum(D(:,j).^2);%求列和
        B(i,j)=-0.5*(D(i,j)^2-dist_td-dist_dj+dist_tdd);%利用MDS算法计算B
    end
end
[V,Lamda]=eig(B)
[LamdaPX,I]=sort(diag(Lamda),'descend')

r=1;
pa=LamdaPX(r);
while sum(pa)/sum(LamdaPX)<1-10^(-15)     %根据特征值所占的比选择出r
      r=r+1;
      pa=[pa LamdaPX(r)];
end
Vs=V(:,I(1:r));
Lamdas=Lamda(I(1:r),I(1:r));
Z=(Lamdas.^1/2)*Vs';