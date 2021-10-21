%% sfls.m

%% Compute the output(s) of an interval singleton type-2 FLS when the 
%% antecedent membership functions are Gaussian primary membership 
%% functions with uncertain means.

%% M1,M2,sigma are mxn matrix denotes the mean and std of
%% antecedent Gaussian MFs (m rules, with n antecedent in each rule)
%% c1,c2 are mx1 vectors, which denotes the COS of consequents
%% X is input matrix, Lxn matrix, each row is an input.
%% R1,R2, and R are the lower, upper and average output

%% call the function gausstype2.m and interval_wtdavg.m


function [R1,R2,R]=sfls_type2(X,M1,M2,sigma,c1,c2);
R1=[];
R2=[];
[L,n]=size(X);
[m,n]=size(M1);
c0=(c1+c2)/2;
s=c2-c0;
%parpool(6)
parfor i=1:L
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
h=(UU+LL)/2;
delta=UU-h;
[l_out,r_out] = interval_wtdavg(c0',s',h,delta);
R1=[R1,l_out];
R2=[R2,r_out];
end

R=(R1+R2)/2;