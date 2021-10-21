close all
clear
clc

x=0:0.001:10;
m1 = 4.2;
m2=6;
s=1;
for i=1:length(x)
[y1(i),y2(i)]=compute_gauss2(x(i),m1,m2,s);
end
plot(x,[y1;y2])
fill([x,x],[y1,y2],'-')
axis([0,10,0,1.2])
hold on
x1=[m1,m1];
[y11,y22]=compute_gauss2(m1,m1,m2,s);
y1=[0,y11];
x2=[m2,m2];
[y111,y222]=compute_gauss2(m2,m1,m2,s);
[y41,y42]=compute_gauss2((m1+m2)/2,m1,m2,s);
x33=[(m1+m2)/2,(m1+m2)/2];
y33=[0,y42];
y2=[0,y111];
plot(x1,y1,'k--')
plot(x2,y2,'k--')
plot(x33,y33,'k--')
plot([0,10],[1,1],'k--')
xlabel('x')
ylabel('\mu(x)')
title('高斯型主隶属函数')





