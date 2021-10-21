clear,clear all
close all
clc
x=load('ex3x.dat');
y=load('ex3y.dat');
x=[ones(length(x),1) x];
n=length(x(1,:));
sigma=std(x);
mu=mean(x);
x(:,2)=(x(:,2)-mu(2))./sigma(2);
x(:,3)=(x(:,3)-mu(3))./sigma(3);

alpha=[0.01,0.03,0.1,0.3,1,1.3];
m=length(y);
T=100;
THETA=[];

for j=1:length(alpha)
    theta=zeros(n,1);
   for i=1:T
       h=x*theta;
       J(i)=(1/(2*m)).*(h-y)'*(h-y);
     %  a=alpha(j);
       theta=theta-alpha(j).*(1/m)*x'*(h-y);
   end
   plot(0:99,J(1:100))
   hold on
end
legend('0.01','0.03','0.1','0.3','1','1.3')
xlabel('iteration number')
ylabel('J')

figure
alpha=1.4;
m=length(y);
T=100;
THETA=[];
theta=zeros(n,1);
for i=1:T
    h=x*theta;
    J(i)=(1/(2*m)).*(h-y)'*(h-y);
 %  a=alpha(j);
    theta=theta-alpha.*(1/m)*x'*(h-y);
end
plot(0:99,J(1:100))
legend('1.4')
xlabel('iteration number')
ylabel('J')

    
    
    
    
    
    