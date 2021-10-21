q1=length(y);
q2=length(g);
X=[y(2:q1-2) y(1:q1-3) g(2:q2-1)'];     %删除规则时所用的数据集
Y=y(3:q1-1);
[m2,n2]=size(X);
for j=1:m2    %计算P
    zz=[];
    for k=1:M       
        t3=[];
        for i=1:n2    
            tt3=exp(-(X(j,i)-mu(k,i))^2/sigma(k,i))^2;
            t3=[t3 tt3];
        end
        zz(k)=prod(t3); 
    end
    b(j)=sum(zz);
    for k=1:M  
        P(j,k)=zz(k)./b(j);
    end
end