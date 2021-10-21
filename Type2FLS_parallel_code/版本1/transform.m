function [M1,M2,sigma,c1,c2] = transform(input,num1,num2)
A=reshape(input(1:(3*num1*num2)),num1,[]);
A=A';
M1=A(1:num2,:);
M2=A(num2+1:2*num2,:);
sigma=A(2*num2+1:3*num2,:);
c1=input(3*num2*num1+1:3*num2*num1+num2)';
c2=input(3*num2*num1+num2+1:3*num2*num1+2*num2)';
end

