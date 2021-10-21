function [ output_x,output_value ] = change_variable_method( c1,c0,d1,d0,A,b )
%change_variable_method(变量变换法)用于求解max{F(x)=(c1*x+c0)/(d1*x+d0)|Ax<=b,x>=0}
[m,n] = size(A);
%% 首先判断分母(d1*x+d0)的符号
[x1,fval1,exitflag1,output1,lambda1] = linprog(d1,A,b,[],[],zeros(length(d1),1));
denominator_min = fval1+d0;
[x2,fval2,exitflag2,output2,lambda2] = linprog(-d1,A,b,[],[],zeros(length(d1),1));
denominator_max=fval2+d0;
%% 根据denominator_min和denominator_max的符号选择方法
c11 = [c1,c0];
c22 = [c1,c0];
A(:,n+1) = -b;
Aeq = [d1,d0];  
beq=1;
if denominator_min >= 0      
    [x11,fval11,exitflag11,output11,lambda11] = linprog(-c11,A,zeros(m,1),Aeq,beq,zeros(length(c11),1));
    output_x = x11(1:length(c1))./x11(length(c11));
    output_value = -fval11;
elseif denominator_max <=0
    [x22,fval22,exitflag22,output22,lambda22] = linprog(c22,A,zeros(m,1),-Aeq,beq,zeros(length(c22),1));
    output_x = x22(1:length(c11)-1)./x22(length(c11));
    output_value = -fval22; 
else
     [x11,fval11,exitflag11,output11,lambda11] = linprog(-c11,A,zeros(m,1),Aeq,beq,zeros(length(c11),1));
     [x22,fval22,exitflag22,output22,lambda22] = linprog(c22,A,zeros(m,1),-Aeq,beq,zeros(length(c22),1));
     x33 = [x11 x22];
     fval33 = [fval11,fval22];
     la1 = find(fval33 == max(fval33));
     output_value = max(fval33);
     outx = x33(la1(1));
     output_x = outx(1:length(c11)-1)./outx(length(c11));
end

end

