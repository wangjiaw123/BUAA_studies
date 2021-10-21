clear,close all
clc
x=load('ex5Linx.dat');
y=load('ex5Liny.dat');
figure
plot(x,y,'ko')
xlabel('x')
ylabel('y')
m=length(x);
n=6;
x=[ones(m,1) x x.^2 x.^3 x.^4 x.^5];
lamda=[0,1,10];
theta=zeros(n,1);
e=eye(n);
e(1,1)=0;
for i=1:length(lamda)
    theta=(x'*x+lamda(i)*e)\x'*y;
    theta_norm=norm(theta)
    xx=linspace(min(x(:,2))-0.1,max(x(:,2))+0.1,50);
    xx=xx';
    xx=[ones(length(xx),1) xx xx.^2 xx.^3 xx.^4 xx.^5];
    yy=xx*theta;
    figure
    plot(x(:,2),y,'ko')
    xlabel('x')
    ylabel('y')
    hold on
    plot(xx(:,2),yy);
    
end


    