function [r_out,I2l,I2u,wr]= rightpoint(c2,LL,UU);

N=length(c2);
[b2,I2]=sort(c2);
L=[];
U=[];

for i=1:N
l=LL(I2(i));
r=UU(I2(i));
L=[L,l];
U=[U,r];
end

r_out=b2*U'/sum(U);

s=b2*U';
s1=sum(U);

I2l=[];
I2u=[];
wr=[];

for i=1:N
s=s-b2(i)*U(i)+b2(i)*L(i);
s1=s1-U(i)+L(i);
 if s/s1>r_out
   r_out=s/s1;
%  wr=[wr,L(i)];
I2l=[I2l,I2(i)];
else
I2u=[I2u,I2(i)];
% wr=[wr,U(i)];
end
end

L=length(I2l);
for i=1:L
wr(I2l(i))=LL(I2l(i));
end

L=length(I2u);
for i=1:L
wr(I2u(i))=UU(I2u(i));
end
