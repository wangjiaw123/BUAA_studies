clear,clear all
close all
clc
x=load('ex2x.dat');
y=load('ex2y.dat');
figure
plot(x,y,'.')
xlabel('age')
ylabel('height')

alpha=0.07;
m=length(y);
x=[ones(length(x),1) x];
T=2000;
theta=[0 0]';
for i=1:T
    h=x*theta;
    theta=theta-alpha/m*x'*(h-y);
end
theta1=linspace(-10,10,100);
theta2=linspace(-10,10,100);
J=zeros(100);
for i=1:length(theta1)
    for j=1:length(theta2)
        t=[theta1(i) theta2(j)];
        J(i,j)=1/(2*m)*(t*theta-y)'*(t*theta-y);
    end
end
figure
surf(theta1,theta2,J)

figure
plot(x,y,'.')
xlabel('age')
ylabel('height')
hold on
xx=2:0.01:8;
xx=xx';
xx=[ones(length(xx),1) xx];
yy=xx*theta;
plot(xx,yy,'k')
axis([2 8 0.7 1.4]);

figure
contour(theta1,theta2,J,logspace(-2,2,15))




