function [R1,R2,R]=sfls_type2_parfor(X,input,M,n)

[ M1,M2,sigma,c1,c2 ] = Mamdani_row_to_discrete( input,M,n );
sigma=sigma+0.001;
% c1=c1+0.001;
% c2=c2+0.001;
R1=[];
R2=[];
[L,n]=size(X);
[m,n]=size(M1);

c0=(c1+c2)/2;
s=c2-c0;

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