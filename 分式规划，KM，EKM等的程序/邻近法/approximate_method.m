function [ fvalue ] = approximate_method( c,d,A,b )
%使用approximate_method(邻近算法)求解max{F(x)=(c*x)/(d*x)|Ax=b,x>=0}
[sizeA_r, sizeA_c] = size(A);
if sizeA_r ~= rank(A)
    disp('矩阵A行不满秩，方程出现冗余!')
    return;
end
%% 找到x的一个可行解
findbaseA = rref(A);                       %将矩阵化成行最简形
cow_location = [];                         %cow_location是矩阵A的极大线性无关组的列位置
[lenfindbaseA_r,lenfindbaseA_c] = size(findbaseA);
for i = 1:lenfindbaseA_c
    findlocction = find(findbaseA(:,i) == 1);
    if (length(findlocction) == 1) && isempty(cow_location)
       cow_location = [i];
    elseif (length(findlocction) == 1) && (length(cow_location) ~= 0)
        gzB = zeros(lenfindbaseA_r,length(cow_location));
        gzB(findlocction,:) = ones(1,length(cow_location));
        if rank(findbaseA(:,cow_location)-gzB) == rank(findbaseA(:,cow_location))
            cow_location = [cow_location i];
        end    
    end
end
xx_1 = A(:,cow_location)\b;
X=zeros(sizeA_c,1);
for i = 1:length(cow_location)
    X(cow_location(i)) = xx_1(i);  %找到初始可行点
end
%% 
lambda1 = (c*X)/(d*X);
[X1star,fval1,exitflag1,outpt1,lama1] = linprog(-d,A,b,c-lambda1*d,0,zeros(length(c),1));
if d*X == d*X1star
    [X1star,fval2,exitflag2,outpt2,lama2] = linprog(d,A,b,c-lambda1*d,0,zeros(length(c),1));
end
X2star=(X1star+X)/2;
%% 
[XX,fval3,exitflag3,outpt3,lama3] = linprog(-c,A,b,lambda1*c+d,(lambda1*c+d)*X2star,zeros(length(c),1));
lambda2=(c*XX)/(d*XX);
while ~(abs(lambda1-lambda2) <= 10e-6)
    X = XX;
    lambda1 = (c*X)/(d*X);
    [X1star,fval1,exitflag1,outpt1,lama1] = linprog(-d,A,b,c-lambda1*d,0,zeros(length(c),1));
    if d*X == d*X1star
       [X1star,fval2,exitflag2,outpt2,lama2] = linprog(d,A,b,c-lambda1*d,0,zeros(length(c),1));
    end
    X2star=(X1star+X)/2;
    [XX,fvalue111,exitflag3,outpt3,lama3] = linprog(-c,A,b,lambda1*c+d,(lambda1*c+d)*X2star,zeros(length(c),1));
    lambda2=(c*XX)/(d*XX); 
end
fvalue=(c*XX)/(d*XX);







