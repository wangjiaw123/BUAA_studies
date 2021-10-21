%��Frank-Wolfe�㷨����ģ��ϵͳ��ʹ���ʶ��ַ���y(k+1)=0.3*y(k)+0.6*y(k-1)+g[u(k)],
%����g(u)=0.6sin(pi*u)+0.3sin(3*pi*u)+0.1sin(5*pi*u),u(k)=cos(2*pi*k/25)
%ģ��ϵͳ��3����1����ģ�ʹ��

clear
close all
clc

%% ��������
y=rand(2,1);
u(1)=sin(2*pi/200);
g(1)=0.6*sin(pi*u(1))+0.3*sin(3*pi*u(1))+0.1*sin(5*pi*u(1));
for k=2:1002
    u(k)=sin(2*pi*k/200);
    g(k)=0.6*sin(pi*u(k))+0.3*sin(3*pi*u(k))+0.1*sin(5*pi*u(k));
    y(k+1)=0.3*y(k)+0.6*y(k-1)+g(k);
end
plot(y)
legend('ʵ��ͼ��')

M=10;   %����ģ��ϵͳ���������ģ��������3����1����ģ�
I=3;
J=1;
% xtrain=[y(2:201) y(1:200) g(2:201)'];   %ѵ��ģ��ϵͳ����������
% ytrain=y(3:202);

xtrain=[y(2:601) y(1:600) g(2:601)'];   %ѵ��ģ��ϵͳ����������
ytrain=y(3:602);
%% ���漸����ʵ��Frank-Wolfe�㷨
tic%����F-W�㷨���õ�ʱ��
epsilon=10e-6;
x0=10*rand(1,M*(2*I+J));%���ó�ʼ�� 
x0_save=x0;%�����ʼ��
Tmax=100;%��������������

rho=0.5;%%%%%%%
alpha=0.4;%%%%% ������Armijo�������������Ĳ���


%% ����������Թ滮ʱ���������½�
% c_mu1=-5*ones(1,40);
% c_mu2=5*ones(1,40);
% for i=1:30
%     if mod(i,3)==0;
%        c_mu1(i+10)=-1; 
%        c_mu2(i+10)=1;
%     end
% end
% sigma_ini1=0.1*ones(1,30);
% sigma_ini2=2*ones(1,30);
% bound_low=[c_mu1 sigma_ini1]';
% bound_up=[c_mu2 sigma_ini2]';

c_mu1=-5*ones(1,40);
c_mu2=5*ones(1,40);
for i=1:30
    if mod(i,3)==0;
       c_mu1(i+10)=-2; 
       c_mu2(i+10)=2;
    end
end
sigma_ini1=0*ones(1,30);
sigma_ini2=2*ones(1,30);
bound_low=[c_mu1 sigma_ini1]';
bound_up=[c_mu2 sigma_ini2]';
%% 
Error=[];
yk_save=[];
t=0;
while t<Tmax %���Ƶ�������Ϊ**
    [ df, Err,mu_df,c_df,sigma_df] = compute_df_f( xtrain,ytrain,x0,M,I,J );
    Error=[Error Err];
    yk=linprog(-df,[],[],[],[],bound_low,bound_up);
    yk_save=[yk_save yk];
    d=yk-x0';
    E1=abs(df*d);
    if(E1<epsilon)
         break; 
    end
%      if(Err<epsilon)
%           break; 
%      end   
    
    m=0;
    while m< 20 %�����Թ滮��Armijo����
        [df,Err]=compute_df_f( xtrain,ytrain,x0,M,I,J );
        fvalue0=Err;
        x0=x0+rho^m*d'; 
        [df1,Err]=compute_df_f( xtrain,ytrain,x0,M,I,J );
        fvalue1=Err;
        if fvalue1<fvalue0+alpha*rho^m*df*d
             mk=m; 
             break;
        end
    end
    x0=x0+rho^mk*d';
    t=t+1;
end
toc
figure
plot(Error)
legend('��ʧ(Ŀ��)����')

x_predict=[y(2:1001) y(1:1000) g(2:1001)'];   %ѵ��ģ��ϵͳ����������
y_predict=y(3:1002);
[ df, Err,fls_f,mu_df,c_df,sigma_df] = compute_df_f( x_predict,...
                      y_predict,x0,M,I,J );%�˴�����compute_df_f������ԭ��
figure                                     %    compute_df_f�������в����Ǽ���ģ��ϵͳ�����fls_f����Ҫʹ��fls_f
plot(fls_f)
hold on
plot(y)
legend('Ԥ��','ʵ��ͼ��')
title('��F-W�㷨�Ż�ģ��ϵͳ�������Ԥ�������ʵ��ͼ��')
figure
plot(fls_f)
legend('Ԥ��')
