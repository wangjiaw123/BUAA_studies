clc,clear
close
clear all
tic
%生成数据集
y=rand(2,1);
u(1)=sin(2*pi/200);
g(1)=0.6*sin(pi*u(1))+0.3*sin(3*pi*u(1))+0.1*sin(5*pi*u(1));
for k=2:602
    u(k)=sin(2*pi*k/200);
    g(k)=0.6*sin(pi*u(k))+0.3*sin(3*pi*u(k))+0.1*sin(5*pi*u(k));
    y(k+1)=0.3*y(k)+0.6*y(k-1)+g(k);
end

M=30;   %假设生成M条规则
xtrain=[y(2:401) y(1:400) g(2:401)'];   %训练模糊系统的数据
[m,n]=size(xtrain);
ytrain=y(3:402);         
sigma=load('sigma.mat'); %初始化高斯型隶属函数的宽度%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma=sigma.sigma(:,1:3);
mu=load('mu.mat');    %初始化高斯型隶属函数的中心
mu=mu.mu(:,1:3);                               %先随机设置初始数据，当得到比较好的结果后得到sigma,mu,yzx,再将它们作为初始数据
yzx=load('yzx.mat');   %初始化后件中心
yzx=yzx.yzx(:,1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sigma=rand(M,n); %初始化高斯型隶属函数的宽度
% mu=rand(M,n);    %初始化高斯型隶属函数的中心
% yzx=rand(M,1);   %初始化后件中心
T=200;
Err=0.05;
count=0;
z=[];
alpha=0.5;
f=zeros(m,1);
while count < T      %用前400组数据和梯度下降法训练参数sigma,mu,yzx
   for j=1:m 
%        for k=1:M        %计算f
%            t1=[];
%            for i=1:n    
%                tt1=exp(-(xtrain(j,i)-mu(k,i))^2/sigma(k,i))^2;
%                t1=[t1 tt1];
%            end
%           z(k)=prod(t1); 
%           a(k)=yzx(k)*z(k); 
%        end
%        b=sum(z);
%        a=sum(a);
%         f(j)=a/b;
         [ff,b,z]  = computer_f( xtrain(j),mu,sigma,yzx,n,M ) ;
         f(j)=ff;


       for k=1:M   %用梯度下降法更新参数
           yzx(k)=yzx(k)-alpha*(f(j)-ytrain(j))/b*z(k);
           for i=1:n
               mu(k,i)=mu(k,i)-alpha*(f(j)-ytrain(j))/b*(yzx(k)-f(j))*z(k)*2*(xtrain(j,i)-mu(k,i))/sigma(k,i)^2;
               sigma(k,i)=sigma(k,i)-alpha*(f(j)-ytrain(j))/b*(yzx(k)-f(j))*z(k)*2*((xtrain(j,i)-mu(k,i))^2)/sigma(k,i)^3;
           end
       end  
   end
   count=count+1;  
%    for k=1:M        %计算f
%        t1=[];
%        for i=1:n    
%            tt1=exp(-(xtrain(j,i)-mu(k,i))^2/sigma(k,i))^2;
%            t1=[t1 tt1];
%        end
%            z(k)=prod(t1); 
%            a(k)=yzx(k)*z(k); 
%    end
%    b=sum(z);
%    a=sum(a);
%    f(j)=a/b;
    [ff,b,z]  = computer_f( xtrain(j),mu,sigma,yzx,n,M ) ;
    f(j)=ff;
   E=0.5*(f(j)-y(j))^2;
   if E<Err
      break;
   end
end

xtest=[y(403:602) y(402:601) g(403:602)'];%生成测试集
[m1,n1]=size(xtest);
   for j=1:m1 
%        for k=1:M        %计算f
%            t1=[];
%            for i=1:n    
%                tt1=exp(-(xtest(j,i)-mu(k,i))^2/sigma(k,i))^2;
%                t1=[t1 tt1];
%            end
%           z(k)=prod(t1); 
%           a(k)=yzx(k)*z(k); 
%        end
%        b=sum(z);
%        a=sum(a);
%        f(j+m)=a/b;
    [ff,b,z]  = computer_f( xtest(j),mu,sigma,yzx,n,M ) ;
    f(j+m)=ff;
   end
for i=1:600
    error(i)=f(i)-y(i);
end
figure
subplot(2,1,1)
plot(f)
hold on
plot(y)
legend('仿真','原始')
subplot(2,1,2)
plot(error)
legend('误差')
e_BP=error;
value_BP=f(1:600)';
