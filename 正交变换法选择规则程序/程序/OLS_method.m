compute_P;
%r=25;%������Ҫ���µĹ�����Ŀ(�����޸�)
U=[];
u=[];
G=[];
gg=[];
I_OLS=[];
err=[];
for i=1:M        %��һ������u(1)��g(1),������Ӧ���д�������I
    U(:,i)=P(:,i);
    G(i)=((U(:,i))'*Y)/((U(:,i))'*U(:,i));
    err(i)=(((G(i))^2)*U(:,i)'*U(:,i))/(y'*y);   
end
location=min(find(err==max(err)));
I_OLS=[I_OLS location];
u(:,1)=U(:,location);
gg(1)=G(location);

II=[1:M];
[r1,c1]=size(II);
jhx=II(r1,1:c1-1);        
for i=location:c1-1
    jhx(:,i)=II(:,i+1);%�ҵ�i1�󣬴�II��ɾȥ��Ӧ���б�
end
II=[];
II=jhx;
alpha=zeros(r);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALPHA=zeros(r);
AL=zeros(r);
for k=2:r
    [r2,c2]=size(II);
    U=[];
    UUU=[];
    G=[];
    GG=[];
    err=[];
    for h=1:c2
        for j=1:k-1
            ALPHA(j,k)=(u(:,j)'*P(:,II(h)))/(u(:,j)'*u(:,j));
        end
        aa=zeros(m2,1);
        for s=1:k-1
            aa=aa+ALPHA(s,k).*u(:,s);  
        end
        U(:,h)=P(:,II(h))-aa;
        G(:,h)=(U(:,h)'*Y)/(U(:,h)'*U(:,h));
        err(h)=((G(:,h)^2)*U(:,h)'*U(:,h))/(Y'*Y);
        AL(:,h)=ALPHA(:,k); 
        GG=[GG G(:,h)];
        UUU=[UUU U(:,h)];
        err=[err err(h)];
    end
    location=min(find(err==max(err)));%���ظ�ȡ��С
    alpha(:,k)=AL(:,location);
    u=[u UUU(:,location)];
    gg=[gg GG(:,location)];
    I_OLS=[I_OLS II(location)];
    [r1,c1]=size(II);
    jhx=II(r1,1:c1-1);        
    for i=location:c1-1
        jhx(:,i)=II(:,i+1);%�ҵ�i1�󣬴�II��ɾȥ��Ӧ���б�
    end
    II=[];
    II=jhx;
end
alpha=eye(r)+alpha;
theta1=alpha\gg';
    
SIGMA1=sigma(I_OLS,:);
MU1=mu(I_OLS,:);

for j=1:m2
    for k=1:r        %����f
        t1=[];
        for i=1:n2    
            tt1=exp(-(X(j,i)-MU1(k,i))^2/SIGMA1(k,i))^2;
            t1=[t1 tt1];
        end
        z(k)=prod(t1); 
        a(k)=theta1(k)*z(k); 
     end
     b=sum(z);
     a=sum(a);
     ff1(j)=a/b;
end

for i=1:600
    error1(i)=ff1(i)-y(i);
end

%figure
subplot(2,1,1)
plot(ff1)
hold on
plot(y)
legend('OLS','ԭʼ')
subplot(2,1,2)
hold on
plot(error1)
legend('Err OLS')
e_OLS=error1;
value_OLS=ff1(1:600)';

























