function [ ff,b,z ] = computer_f( x,mu,sigma,yzx,n,M )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
a=0;b=0;
for k=1:M        %计算f
    t1=[];
    for i=1:n    
        tt1=exp(-(x-mu(k,i))^2/sigma(k,i))^2;
        t1=[t1 tt1];
    end
        z(k)=prod(t1); 
        a(k)=yzx(k)*z(k); 
end
b=sum(z);
a=sum(a);
ff=a/b;
end

