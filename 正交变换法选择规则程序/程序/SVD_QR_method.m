compute_P;
%r=25;%设置想要留下的规则数目(可以修改)
ysq=y(1:q1-3);
[Usvd,Ssvd,Vsvd]=svd(P);
Aqr=[Vsvd(1:r,1:r)', Vsvd((r+1):M,1:r)'];
[Q,R,E]=qr(Aqr);
PP=P*E;
SQ=[];
for i=1:M
   location=find(E(:,i)==1);
   SQ=[SQ location];
end
I_SVD_QR=SQ(1:r);
Pr=PP(:,1:r);
theta3=(Pr'*Pr)\(Pr'*ysq);

SIGMA3=sigma(I_SVD_QR,:);
MU3=mu(I_SVD_QR,:);

for j=1:m2
    for k=1:r        %计算f
        t1=[];
        for i=1:n2    
            tt1=exp(-(X(j,i)-MU3(k,i))^2/SIGMA3(k,i))^2;
            t1=[t1 tt1];
        end
        z(k)=prod(t1); 
        a(k)=theta3(k)*z(k); 
     end
     b=sum(z);
     a=sum(a);
     ff3(j)=a/b;
end

for i=1:600
    error1(i)=ff3(i)-y(i);
end

% figure
subplot(2,1,1)
plot(ff3)
hold on
plot(y)
legend('SVD QR','原始')
subplot(2,1,2)
hold on
plot(error1)
legend('Err SVD QR')
e_SVD_QR=error1;
value_SVD_QR=ff3(1:600)';





