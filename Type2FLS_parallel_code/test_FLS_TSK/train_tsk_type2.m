%% train_tsk_type2.m

%% Tune the parameters of an interval type-2 TSK FLS A2-C1 when the 
%% antecedent membership functions are Gaussian primary membership 
%% functions with uncertain means, using some input–output training data.

%% M1,M2, sigma are mxn matrix denotes the mean and std of
%% antecedent Gaussian MFs (N rules, with n antecedent in each rule)
%% C, S are Nx(ant+1) matrix denoting the center and spread of consequents para.
%% X is input matrix, L(x)ant matrix, each row is onw input.
%% D is Lx1 vector which denotes the desired output


function [M1,M2,C,S,sigma]=train_sfls_type2(x,D,M1,M2,sigma,C,S,alpha);
X=x;
[L,n]=size(x);
[N,n]=size(M1);

for i=1:L
%U=[];
%MU1=[];
UU=[];
LL=[];
for j=1:N
Uu=1;
Ll=1;
for m=1:n
P=[sigma(j,m),M1(j,m),M2(j,m)];

[uu,ll]=gausstype2(x(i,m),P);
Uu=Uu*uu;
Ll=Ll*ll;
end
UU=[UU,Uu];
LL=[LL,Ll];
end

c2=C(:,1)+S(:,1);
c1=C(:,1)-S(:,1);
for t=1:length(x(1,:))
c2=c2+C(:,t+1)*x(i,t)+S(:,2)*abs(x(i));
c1=c1+C(:,t+1)*x(i,t)-S(:,2)*abs(x(i));
end

[r_out,I2l,I2u,wr]= rightpoint(c2',LL,UU);
[l_out,I1u,I1l,wl]= leftpoint(c1',LL,UU);

f=(l_out+r_out)/2;
e=D(i)-f;


fa1=wl/sum(wl);
fa2=wr/sum(wr);

%%% Select MFs contributed to the right-most 

ME10=M1;
ME20=M2;
sigma0=sigma;

LE=length(I2u);

for t=1:LE
for k=1:n
if x(i,k)<M1(I2u(t),k)
l=I2u(t);
M1(l,k)=M1(l,k)+alpha*e*0.5*((x(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(I2u(t))/sum(wr);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((x(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
elseif  x(i,k)>M2(I2u(t),k)
l=I2u(t);
M2(l,k)=M2(l,k)+alpha*e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(l)/sum(wr);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
end
end
end


LE=length(I2l);

for t=1:LE
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

%%% Select MFs contributed to the left-most 

LE=length(I1l);

for t=1:LE
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

LE=length(I1u);

for t=1:LE
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

sigma=sigma0;

for l=1:N
for k=1:n
if sigma(l,k) < 0
sigma(l,k)=abs(sigma(l,k));
end
end
end

fa1=wr'/sum(wr);
fa2=wl'/sum(wl);


C(:,1)=C(:,1)+alpha*e*(fa1+fa2)/2;
S(:,1)=S(:,1)+alpha*e*(fa2-fa1)/2;
for j=2:(n+1)
C(:,j)=C(:,j)+alpha*e*X(i,j-1)*(fa1+fa2)/2;
S(:,j)=S(:,j)+alpha*e*abs(X(i,j-1))*(fa2-fa1)/2;
end

%%%% After Training, The type-2 MFs should be reconstructed

S=abs(S);

ME1=[];
ME2=[];
for t=1:N
P=[M2(t,:)',M1(t,:)']';
m2=max(P);
m1=min(P);
ME1=[ME1,m1'];
ME2=[ME2,m2'];
end
M1=ME1';
M2=ME2';

end  




end