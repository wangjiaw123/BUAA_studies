compute_P;
[ U,UP,S,SP,V,VP,I,E,count ] = onesideJacobi1( P );
%r=25;%设置想要留下的规则数目(可以修改)
II=I(1:r);
Pr=P(:,II);
theta5=(Pr'*Pr)\Pr'*Y;
SIGMA5=sigma(II,:);
MU5=mu(II,:);

for j=1:m2
    for k=1:r        %计算f
        t1=[];
        for i=1:n2    
            tt1=exp(-(X(j,i)-MU5(k,i))^2/SIGMA5(k,i))^2;
            t1=[t1 tt1];
        end
        z(k)=prod(t1); 
        a(k)=theta5(k)*z(k); 
     end
     b=sum(z);
     a=sum(a);
     ff5(j)=a/b;
end

for i=1:600
    error1(i)=ff5(i)-y(i);
end

% figure
subplot(2,1,1)
plot(ff5)
hold on
plot(y)
legend('DirectSVD','原始')
subplot(2,1,2)
hold on
plot(error1)
legend('Err DSVD')
e_DSVD=error1;
value_DSVD=ff5(1:600)';

