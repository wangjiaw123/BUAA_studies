compute_P;
%r=25;%设置想要留下的规则数目(可以修改)
Y=y(1:q1-3);
Py=[P,Y];
[UW,SW,VW]=svd(Py);
Aqr4=[VW(1:r,1:r)', VW((r+1):M,1:r)'];
[QE,RE,EW]=qr(Aqr4);
PPW=P*EW;
SQW=[];
for i=1:M
   location=find(EW(:,i)==1);
   SQW=[SQW location];
end
I_TLS=SQW(1:r);
Prw=PPW(:,1:r);
[UJ,SJ,VJ]=svd([Prw Y]);
theta4=(-1/VJ(r+1,r+1)).*VJ(:,r+1);

SIGMA4=sigma(I_TLS,:);
MU4=mu(I_TLS,:);

for j=1:m2
    for k=1:r        %计算f
        t1=[];
        for i=1:n2    
            tt1=exp(-(X(j,i)-MU4(k,i))^2/SIGMA4(k,i))^2;
            t1=[t1 tt1];
        end
        z(k)=prod(t1); 
        a(k)=theta4(k)*z(k); 
     end
     b=sum(z);
     a=sum(a);
     ff4(j)=a/b;
end

for i=1:600
    error1(i)=ff4(i)-y(i);
end

% figure
subplot(2,1,1)
plot(ff4)
hold on
plot(y)
legend('TLS','原始')
subplot(2,1,2)
hold on
plot(error1)
legend('Err TLS')
e_TLS=error1;
value_TLS=ff4(1:600)';