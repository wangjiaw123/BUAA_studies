function [ M1,M2,sigma,c1,c2 ] = Mamdani_row_to_discrete( input,M,n )
%row_to_discrete 将行向量(参数混合)转换为分开的参数矩阵
input=reshape(input,[1,3*M*n+2*M]);
M1=reshape(input(1:M*n),[n,M])';
M2=reshape(input(1+M*n:2*M*n),[n,M])';
sigma=reshape(input(1+2*M*n:3*M*n),[n,M])';
c1=reshape(input(1+3*M*n:3*M*n+M),[M,1]);
c2=reshape(input(M+1+3*M*n:3*M*n+2*M),[M,1]);
end

