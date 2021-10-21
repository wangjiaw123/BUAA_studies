clear,clc
m=10;%����ԭʼ�ռ��ά��
r=5;%���õ�ά�ռ��ά��
D=10*rand(m);%�������һ������
a=0;
B=zeros(m);
for i=1:m
    for j=1:m
        a=a+D(i,j)^2;
    end
end
dist_tdd=(1/m^2)*a;
for i=1:m
    for j=1:m
        dist_td=(1/m)*sum(D(i,:).^2);
        dist_dj=(1/m)*sum(D(:,j).^2);
        B(i,j)=-0.5*(D(i,j)^2-dist_td-dist_dj+dist_tdd);
    end
end
[V,Lamda]=eig(B)
[LamdaPX,I]=sort(diag(Lamda),'descend')
Vs=V(:,I(1:r));
Lamdas=Lamda(I(1:r),I(1:r));
Z=(Lamdas.^1/2)*Vs';
