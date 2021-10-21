function [M1,M2,c1,c2,sigma,I2l,I2u,I1u,I1l]=train_sfls_type2(X,D,M1,M2,sigma,c1,c2,alpha)

[L,n]=size(X);
[N,n]=size(M1);
M12=M1;M22=M2;sigma22=sigma;
for i=1%:L
U=[];
MU1=[];
UU=[];
LL=[];
for j=1:N
Uu=1;
Ll=1;
for m=1:n
sigma(j,m)
M1(j,m)
M2(j,m)
P=[sigma(j,m),M1(j,m),M2(j,m)]
i
j
m
[uu,ll]=gausstype2(X(i,m),P);
Uu=Uu*uu;
Ll=Ll*ll;
P=[];
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

LE1=length(I2u);
LE2=length(I2l);
LE3=length(I1l);
LE4=length(I1u);

spmd
    labindex
if labindex==1    
for t=1:LE1
for k=1:n 
if X(i,k)<M1(I2u(t),k)
l=I2u(t);
M1(l,k)=M1(l,k)+alpha*e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(I2u(t))/sum(wr);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
elseif  X(i,k)>M2(I2u(t),k)
l=I2u(t);
M2(l,k)=M2(l,k)+alpha*e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(l)/sum(wr);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
end
end
end
end

% LE2=length(I2l);
if labindex==2 
for t=1:LE2
for k=1:n
if X(i,k)<(M1(I2l(t),k)+M2(I2l(t),k))/2
l=I2l(t);
M2(l,k)=M2(l,k)+alpha*e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(l)/sum(wr);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
else
l=I2l(t);
M1(l,k)=M1(l,k)+alpha*e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(l)/sum(wr);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
end
end
end   
end
%%% Select MFs contributed to the left-most 

% LE3=length(I1l);
if labindex==3
for t=1:LE3
for k=1:n 
if X(i,k)<(M1(I1l(t),k)+M2(I1l(t),k))/2
l=I1l(t);
M2(l,k)=M2(l,k)+alpha*e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
else
l=I1l(t);
M1(l,k)=M1(l,k)+alpha*e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
end
end
end
end
% LE4=length(I1u);
if labindex==4
for t=1:LE4
for k=1:n
if  X(i,k)< M1(I1u(t),k)
l=I1u(t);
M1(l,k)=M1(l,k)+alpha*e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
elseif  X(i,k)> M2(I1u(t),k)
l=I1u(t);
M2(l,k)=M2(l,k)+alpha*e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
end
end
end  
end
end%end spmd
sigma=sigma0;

sigma3=zeros(N,n);
M13=zeros(N,n);
M23=zeros(N,n);
for ii=1:n
   M13=M13+cell2mat(M1(ii));
   M23=M23+cell2mat(M2(ii));
   sigma3=sigma3+cell2mat(sigma(ii));
end
sigma=sigma3-(n-1)*sigma22;
M1=M13-(n-1)*M22;
M2=M23-(n-1)*M12;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for l=1:N
% for k=1:n
% if sigma(l,k) < 0
% sigma(l,k)=abs(sigma(l,k));
% end
% end
% end

%fa1=wr'/sum(wr);
%fa2=wl'/sum(wl);
c1=c1+alpha*e*fa1/2;
c2=c2+alpha*e*fa2/2;      %%%%选用KM方法时用此4行命令，下面4行改成注释
% fa1=wr/sum(wr);
% fa2=wl/sum(wl);
% c1=c1+alpha*e*fa1/2;         %% 选用PSA，变量变换法，迭代法时用此4行命令
% c2=c2+alpha*e*fa2/2;
end
end