clear,close all
clc
x=load('ex5Logx.dat');
y=load('ex5Logy.dat');
postive=find(y==1);
negetive=find(y==0);
figure
plot(x(postive,1),x(postive,2),'r+')
hold on
plot(x(negetive,1),x(negetive,2),'o')
legend('y=1','y=0')
xlabel('x1')
ylabel('x2')

g=inline('1.0./(1+exp(-z))');
xx=map_feature(x(:,1),x(:,2));
[m,n]=size(xx);

T=15;
e=eye(n);
e(1,1)=0;
lamda=[0.1 1 10];
for k=1:length(lamda)
    theta=zeros(n,1);
for i=1:T
    z=xx*theta;
    h=g(z);
    J(i)=(1/m)*sum(-y'*log(h)-(1-y)'*log(1-h))+(lamda(k)/m)*(norm(theta))^2;
    G=(lamda(k)/m).*theta;
    G(1)=0;
    grad=(1/m).*xx'*(h-y)+G;
    L=(lamda(k)/m).*e;
    H=(1/m).*xx'*diag(h)*diag(1-h)*xx+L;
    theta=theta-H\grad; 
end
theta
THETE_norm=norm(theta)
figure
plot(0:T-1,J)
legend('J')
xlabel('iteration number')
ylabel('J')

u=linspace(-1,1.5,100);
v=linspace(-1,1.5,100);
for i=1:length(u)
    for j=1:length(v)
        zz(i,j)=map_feature(u(i),v(j))*theta;
    end
end
zz=zz';
figure
plot(x(postive,1),x(postive,2),'r+')
hold on
plot(x(negetive,1),x(negetive,2),'o')
legend('y=1','y=0')
xlabel('x1')
ylabel('x2')
contour(u,v,zz,[0 0])
end





