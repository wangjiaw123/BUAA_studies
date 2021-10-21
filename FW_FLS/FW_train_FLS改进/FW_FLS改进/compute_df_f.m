function [ df, Err,fls_f,mu_df,c_df,sigma_df] = compute_df_f( xtrain,ytrain,x0,M,I,J )
% compute_df_f ���㺯����һ��ƫ��������ͼ������ֵErr
% xtrain��ѵ�����ݶԵ����벿�֣�ytrain��ѵ�����ݵ��������
% Mָģ��ϵͳ�Ĺ��������I�ǹ���ǰ��(ʹ�ø�˹���������������в�����mu,sigma)
%    �ĸ�����J�ǹ�����c�ĸ�����ʹ������ƽ����ģ����
% x0�ǳ�ʼ��ģ��ϵͳ�Ĳ�����x0�ɿ��������������ʽΪ[c(1),...,c(M),mu(1),...,
%     mu(M*I),sigma(1),...,sigma(M*J)]
%% ��ģ��ϵͳ�Ĳ������и�ֵ
c=x0(1:M)';
mu=reshape(x0(M+1:M+M*I),I,[]);
mu=mu';%
sigma=reshape(x0(M*I+M+1:M*I*2+M),I,[]);
sigma=sigma';
%% ����
len_xtrain=length(xtrain(:,1));
fls_f=zeros(1,len_xtrain);

c_df=zeros(M,J,len_xtrain);
c_df(:,:,1)=c;
mu_df=zeros(M,I,len_xtrain);
mu_df(:,:,1)=mu;
sigma_df=zeros(M,I,len_xtrain);
sigma_df(:,:,1)=sigma;
for k=1:len_xtrain
    z=zeros(1,M);
    a=zeros(1,M);
    for i=1:M
        t1=zeros(1,I);
        for j=1:I
           t1(j)=exp(-((xtrain(k,j)-mu(i,j))/sigma(i,j))^2);
        end  
        z(i)=prod(t1); 
        a(i)=c(i)*z(i);   
    end
    b=sum(z);
    a=sum(a);
    fls_f(k)=a/b;
    
    for i1=1:M   %�ֱ�����ݶȲ�����***_df����
        c_df(i1,:,k)=(fls_f(k)-ytrain(k))*z(i1)/b;
        for i2=1:I
            mu_df(i1,i2,k)=(fls_f(k)-ytrain(k))/b*(c(i1,:)-fls_f(j))*z(i1)*2*(xtrain(k,i2)-mu(i1,i2))/sigma(i1,i2)^2;
            sigma_df(i1,i2,k)=(fls_f(k)-ytrain(k))/b*(c(i1,:)-fls_f(k))*z(i1)*2*((xtrain(k,i2)-mu(i1,i2))^2)/sigma(i1,i2)^3;
        end
    end  
    %k=k+1;
end
Err=0;
c_end=zeros(M,J);
mu_end=zeros(M,I);
sigma_end=zeros(M,I);
for i=1:len_xtrain
   Err=Err+0.5*(fls_f(i)-ytrain(i))^2;
   c_end=c_end+c_df(:,:,i);
   mu_end=mu_end+mu_df(:,:,i);
   sigma_end=sigma_end+sigma_df(:,:,i);
end
x0(1:M)=reshape(c_end,1,[]);
x0(M+1:M*I+M)=reshape(mu_end',1,[]);
x0(M*I+M+1:2*M*I+M)=reshape(sigma_end',1,[]);
df=x0;

end

