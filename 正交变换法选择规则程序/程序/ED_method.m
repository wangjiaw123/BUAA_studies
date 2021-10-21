compute_P;
%r=25;%设置想要留下的规则数目(可以修改)
phi_pp=P'*P;
phi_y=P'*y(1:q1-3);
[V,Lambda]=eig(phi_pp);
I_ED=diag(V(:,1:r)*V(:,1:r)');
[Vb,IED]=sort(I_ED,'descend');
% phi_ppr=phi_pp(IED(1:r),:);
% phi_ppr=phi_pp(IED(1:r),IED(1:r));
% phi_yr=phi_y(IED(1:r),:);
phip=P(:,IED(1:r));
phi_ppr=phip'*phip;
phi_yr=phi_y(IED(1:r),:);

theta2=(phi_ppr'*phi_ppr)\(phi_ppr'*phi_yr);

SIGMA2=sigma(IED(1:r),:);
MU2=mu(IED(1:r),:);

for j=1:m2
    for k=1:r        %计算f
        t1=[];
        for i=1:n2    
            tt1=exp(-(X(j,i)-MU2(k,i))^2/SIGMA2(k,i))^2;
            t1=[t1 tt1];
        end
        z(k)=prod(t1); 
        a(k)=theta2(k)*z(k); 
     end
     b=sum(z);
     a=sum(a);
     ff2(j)=a/(b);
end

for i=1:600
    error1(i)=ff2(i)-y(i);
end

% figure
subplot(2,1,1)
plot(ff2)
hold on
plot(y)
legend('ED','原始')
subplot(2,1,2)
hold on
plot(error1)
legend('Err ED')
e_ED=error1;
value_ED=ff2(1:600)';