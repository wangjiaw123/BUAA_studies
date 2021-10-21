function [ X1,F ] = iterative_method( c1,c0,d1,d0,A,b )
%使用iterative_method(迭代算法)求解max{F(x)=(c1*x+c0)/(d1*x+d0)|Ax=b,x>=0}
z=0;
[X1,fv1,exitflag1,outpt1,lama1] = linprog(-(c1-z*d1),A,b,[],[],zeros(length(c1),1));
z=(c1*X1+c0)/(d1*X1+d0);
Fx=(c1*X1+c0)-z*(d1*X1+d0);
while ~(abs(Fx) <= 10e-4)
    [X1,fv1,exitflag1,outpt1,lama1] = linprog(-(c1-z*d1),A,b,[],[],zeros(length(c1),1));
    z=(c1*X1+c0)/(d1*X1+d0);
    Fx=(c1*X1+c0)-z*(d1*X1+d0);
end
F=z;
end

