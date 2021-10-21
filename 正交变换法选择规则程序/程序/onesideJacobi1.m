function [ U,UP,S,SP,V,VP,I,E,count ] = onesideJacobi1( A )
%onesideJacobi 利用单侧雅克比变换对矩阵A进行SVD分解,并且在令向量两两正交时使用前向扫描的方法
%A：输入的需要进行SVD分解的矩阵
%U,SUM,V指A=U*SUM*V'
%B指B=A*V
%E指任意两列(互不相同)向量做内积的和（与0越接近越好）
%SP是对S排序后得到的
%UP指对U排序后得到的
%VP指对V排序后得到的

[m,n]=size(A);
if m<n
    disp('The number of rows must be larger than the number of columns!')
end
if mod(n,2)==1
    disp('The number of column must be even!')
end

count=0;
Err=1e-20;
J=eye(n);
B=eye(n);
while count<100
        a1=[n-1:-2:1];
        a2=[n:-2:2];
    for step=1:n-1

    for i=1:n/2
        tau(i)=(A(:,a1(i))'*A(:,a1(i))-A(:,a2(i))'*A(:,a2(i)))/(2*A(:,a1(i))'*A(:,a2(i)));
        t(i)=(sign(tau(i)))/(abs(tau(i))+(1+tau(i)^2)^0.5);                   
        c(i)=1/((1+t(i)^2)^0.5);                                  %此处计算n/2个正交矩阵的参数
        s(i)=t(i)*c(i);
    end
    for k=1:n/2
        JJ=eye(n);
        JJ(a1(k),a1(k))=c(k); JJ(a2(k),a2(k))=c(k);JJ(a1(k),a2(k))=-s(k);JJ(a2(k),a1(k))=s(k);
        J=J*JJ; 
        A=A*JJ; 
    end
   B=A;
    if mod(step,2)==1
       location=find(a2==n);
       jhx1=a1(location);                                      %此处用于交换位置
       a1(location)=a2(location);
       a2(location)=jhx1;
       jhx2=a2;
       a2(1)=jhx2(n/2);
       for ii=2:n/2                                           %此处用于错位
           a2(ii)=jhx2(ii-1);
       end    
    end
    if mod(step,2)==0
       location=find(a1==n);
       jhx1=a2(location);                                     %此处用于交换位置
       a2(location)=a1(location);
       a1(location)=jhx1;
       jhx2=a2;
       a2(1)=jhx2(n/2);
       for ii=2:n/2
           a2(ii)=jhx2(ii-1);                               %此处用于错位
       end    
    end 
  end 
  count=count+1;   
  a=[];
   for i=1:n-1
       for j=i+1:n
           b=(B(:,i)'*B(:,j))^2;
           a=[a,b];
       end
   end
 
   E=sum(a);              %计算任意两列(互不相同)向量做内积的和
   if abs(E)<=Err
       break;
   end

end
V=J;                     
U=zeros(m);
for i=1:n
    c=0;
    for j=1:m
        c=c+(B(j,i))^2;         
    end
    c=c^0.5;
    sigma(i)=c;          
   U(:,i)= 1/sigma(i).*(B(:,i));      %对B的每一列进行标准化
end

S=zeros(m,n);
for j=1:n
    S(j,j)=sigma(j);         %得出矩阵S对角线上的元
end
[Skpz,I]=sort(diag(S),'descend');   %对S对角线上的元进行排序，并得到对应的列标向量I
UP=U(:,I);
S2=S(:,I);
for i=1:n                            %对S2进行排序输出SP
    if max(S2(:,i))~=0               %考虑到奇异值有可能为0的情况
       location=find(S2(:,i)~=0);
       if location~=i
          S2(i,i)=S2(location,i);
          S2(location,i)=0;
       end
       if location==i
          S2(i,i)=S2(location,i); 
       end
    end
    if max(S2(:,i))==0
        S2(i,i)=0;
    end
end
SP=S2;
VP=V(:,I);
end