%% tsk_type2.m

%% Compute the output(s) of an interval type-2 TSK FLS A2-C1 (type-2 
%% antecedents and type-1 consequent) when the antecedent membership
%% functions are Gaussian primary membership functions with uncertain
%% means.

%% M1,M2,sigma are mxn matrix denotes the mean and std of
%% antecedent Gaussian MFs (m rules, with n antecedent in each rule)
%% C,S are mx(n+1) matrix, which denotes the center and spread of coefficent of parameters
%% X is input matrix, Lxn matrix, each row is onw input.
%% R1,R2, and R are the lower, upper and average output
%% call the function gausstype2.m and interval_wtdavg.m


function [R1,R2,R]=tsk_type2(X,M1,M2,sigma,C,S);

R1=[];
R2=[];
[L,n]=size(X);
[m,n]=size(M);

c0=C(:,1);
s=S(:,1);
for t=1:n
c0=c0+C(:,t)*X(i,t);
s=s+S(:,t)*abs(X(i,t));
end

%c0=(c1+c2)/2;
%s=c2-c0;

for i=1:L
%U=[];
%MU1=[];

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

h=(UU+LL)/2;
delta=UU-h;

[l_out,r_out] = interval_wtdavg(c0',s',h,delta);

R1=[R1,l_out];
R2=[R2,r_out];
end

R=(R1+R2)/2;
end