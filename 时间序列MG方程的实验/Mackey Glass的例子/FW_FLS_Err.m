function [Err]=FW_FLS_Err(X,D,M1,M2,sigma,c1,c2);
[L,n]=size(X);
[N,n]=size(M1);
M11=zeros(size(M1));
M22=zeros(size(M2));
sigma11=zeros(size(sigma));
c11=zeros(size(c1));
c22=zeros(size(c2));
e=zeros(1,L);
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
e(i)=(D(i)-f)^2/2;
end
Err=sum(e);
end
