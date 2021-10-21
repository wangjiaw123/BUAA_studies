clear,close all
clc
x=load('ex4x.dat');
y=load('ex4y.dat');
postive=find(y==1);
negetive=find(y==0);
figure
plot(x(postive,1),x(postive,2),'+')
hold on
plot(x(negetive,1),x(negetive,2),'*')
legend('postive','negetive')
xlabel('x1')
ylabel('x2')

x=[ones(length(x),1) x];
theta=zeros(length(x(1,:)),1);
m=length(x);
g=inline('1./(1+exp(-z))');
T=8;
J=zeros(T,1);

for  i=1:T
    z=x*theta;
    h=g(z);
    grad=(1/m).*x'*(h-y);
    
    H=(1/m).*x'*diag(h)*diag(1-h)*x;
  %下面是自己写的求海赛矩阵的程序，输出的结果始终有问题，
  %我一直没找出来错在那里，就使用了pdf中所给的求H的代码
  %  n=length(theta);
  %  H=zeros(n);
  %  aa=0;
  %  for t=1:n
  %     for j=1:n
  %         for k=1:m
  %            xx=x(k,:)*theta;
  %            aa=aa+g(xx)*(1-g(xx))*x(k,t)*x(k,j);
  %            H(t,j)=aa/m;  
  %         end
  %     end
  % end
   
    theta=theta-H\grad;
    J(i)=(1/m)*sum(-y.*log(h)-(1-y).*log(1-h));
end
theta
plot_X=[min(x(:,2))-2,max(x(:,2))+2];
plot_Y=(-1./theta(3)).*(theta(2).*plot_X+theta(1));
plot(plot_X,plot_Y)
figure
plot(J)
legend('J')
xlabel('iteration number')
ylabel('J')