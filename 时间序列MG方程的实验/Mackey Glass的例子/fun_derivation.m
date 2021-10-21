function [M11,M22,c11,c22,sigma11]=fun_derivation(X,D,M1,M2,sigma,c1,c2);
[L,n]=size(X);
[N,n]=size(M1);
M11=zeros(size(M1));
M22=zeros(size(M2));
sigma11=zeros(size(sigma));
c11=zeros(size(c1));
c22=zeros(size(c2));
for i=1:L
U=[];
MU1=[];
UU=[];
LL=[];
for j=1:N
Uu=1;
Ll=1;
for m=1:n
P=[sigma(j,m),M1(j,m),M2(j,m)];
[uu,ll]=gausstype2(X(i,m),P);
Uu=Uu*uu;
Ll=Ll*ll;
end
UU=[UU,Uu];
LL=[LL,Ll];
end

%% Compute the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r_out,I2l,I2u,wr]= rightpoint(c2',LL,UU);          %此处用了EKM算法
[l_out,I1u,I1l,wl]= leftpoint(c1',LL,UU);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=(l_out+r_out)/2;
e=D(i)-f;

fa1=wl/sum(wl);
fa2=wr/sum(wr);

%%% Select MFs contributed to the right-most 
ME10=M1;
ME20=M2;
sigma0=sigma;

for it=1:N
c22(it)=0.5*e*wr(it)/sum(wr);
end

for it=1:N
c11(it)=0.5*e*wl(it)/sum(wl);
end

LE=length(I2u);
for t=1:LE
for k=1:n 
if X(i,k)<M1(I2u(t),k)
l=I2u(t);
M11(l,k)=e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(I2u(t))/sum(wr);
sigma0(l,k)=e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
elseif  X(i,k)>M2(I2u(t),k)
l=I2u(t);
M22(l,k)=e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(l)/sum(wr);
sigma0(l,k)=e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
end
% l=I2u(t);
% c22(i)=e*wr(l)/sum(wr);
end
end


LE=length(I2l);

for t=1:LE
for k=1:n
if X(i,k)<(M1(I2l(t),k)+M2(I2l(t),k))/2
l=I2l(t);
M22(l,k)=e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(l)/sum(wr);
sigma0(l,k)=e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
else
l=I2l(t);
M11(l,k)=e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(l)/sum(wr);
sigma0(l,k)=e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
end
% l=I2l(t);
% c11(l)=e*wr(l)/sum(wr);
end
end   

%%% Select MFs contributed to the left-most 

LE=length(I1l);

for t=1:LE
for k=1:n 
if X(i,k)<(M1(I1l(t),k)+M2(I1l(t),k))/2
l=I1l(t);
M22(l,k)=e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
else
l=I1l(t);
M11(l,k)=e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
end
% l=I1l(t);
% c11(l)=e*wr(l)/sum(wr);
end
end

LE=length(I1u);

for t=1:LE
for k=1:n
if  X(i,k)< M1(I1u(t),k)
l=I1u(t);
M11(l,k)=e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
elseif  X(i,k)> M2(I1u(t),k)
l=I1u(t);
M22(l,k)=e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
end
% l=I1u(t);
% c22(l)=e*wr(l)/sum(wr);
end
end   
sigma11=sigma0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for l=1:N
% for k=1:n
% if sigma(l,k) < 0
% sigma(l,k)=abs(sigma(l,k));
% end
% end
% end

% fa1=wr'/sum(wr);
% fa2=wl'/sum(wl);
% c11=e*fa1/2;
% c22=e*fa2/2;     
end