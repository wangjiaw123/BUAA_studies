function [u,l,xx]=gausstype2(x,P)

sigma=abs(P(:,1))+0.000001;
m1=min([P(:,2),P(:,3)]);
m2=max([P(:,3),P(:,2)]);

L=length(x);

u=[];
l=[];

for i=1:L
if x(i)>=m1 && x(i)<=(m1+m2)/2
mu1=1; 
mu2=exp(-0.5*((x(i)-m2)^2)/(sigma^2));

elseif  x(i)>(m1+m2)/2 && x(i)<=m2
mu1=1;
mu2=exp(-0.5*((x(i)-m2)^2)/(sigma^2));
else
mu1=exp(-0.5*((x(i)-m2)^2)/(sigma^2));
mu2=exp(-0.5*((x(i)-m2)^2)/(sigma^2));
end

u=[u,max([mu1,mu2])];
l=[l,min([mu1,mu2])];
end

xx=x;
