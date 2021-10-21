%% svd_qr_sfls.m

%% Rule-reduction of an interval singleton type-2 FLS when the 
%% antecedent membership functions are Gaussian primary membership 
%% functions with uncertain means, using some inputÐoutput training data.

%% M1,M2,sigma are mxn matrix denotes the mean and std of
%% antecedent Gaussian MFs (m rules, with n antecedent in each rule)
%% c1,c2 are mx1 vectors, which denotes the COS of consequents
%% X is input matrix, Lxn matrix, each row is an input.
%% call the function gausstype2.m
%% th is a threshold for choosing the principal components

function [M1,M2,sigma,c1,c2]=sfls_type2(X,M1,M2,sigma,c1,c2);

FBFu=[];
FBFl=[];
[L,n]=size(X);
[m,n]=size(M1);
rule=m;

c0=(c1+c2)/2;
s=c2-c0;

for i=1:L
U=[];
MU1=[];

UU=[];
LL=[];
for j=1:m
Uu=1;
Ll=1;
for t=1:n
P=[sigma(j,t),M1(j,t),M2(j,t)];
[uu,ll]=gausstype2(X(i,t),P);
Uu=Uu*uu;
Ll=Ll*ll;
end
UU=[UU,Uu];
LL=[LL,Ll];
end

[r_out,I2l,I2u,wr]= rightpoint(c2',LL,UU);
[l_out,I1u,I1l,wl]= leftpoint(c1',LL,UU);
FBFu=[FBFu,wr'/sum(wr)];
FBFl=[FBFl,wl'/sum(wl)];
end

FBFu=FBFu';
FBFl=FBFl';

%%%%%%%%%%
[U,SS,V]=svd(FBFu);


a=diag(SS);
for i=1:N
if a(i)>th
r=i;
end
end

A=[V(1:r,1:r)', V((r+1):N,1:r)'];
[Q,R,E]=qr(A);
[e,I]=sort(E);
I1=I(N,1:r);
%%%%%%%%%%%%%
[U,SS,V]=svd(FBFl);

a=diag(SS);
for i=1:N
if a(i)>th
r=i;
end
end

A=[V(1:r,1:r)', V((r+1):N,1:r)'];
[Q,R,E]=qr(A);
[e,I]=sort(E);
I2=I(N,1:r);
%%%%%%%%%%%%
I=union(I1,I2);
r=length(I);

%%%%%%% Keep the most important parameters

M1=M1(I,:);
M2=M2(I,:);
c1=c1(I);
c2=c2(I);
sigma=sigma(I,:);
