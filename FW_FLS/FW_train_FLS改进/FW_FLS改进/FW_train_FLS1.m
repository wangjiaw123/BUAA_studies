%用Frank-Wolfe算法调节模糊系统，使其辨识差分方程y(k+1)=0.3*y(k)+0.6*y(k-1)+g[u(k)],
%其中g(u)=0.6sin(pi*u)+0.3sin(3*pi*u)+0.1sin(5*pi*u),u(k)=cos(2*pi*k/25)
%模糊系统是3输入1输出的，使用

clear
close all
clc

%% 产生数据
y=rand(2,1);
u(1)=sin(2*pi/200);
g(1)=0.6*sin(pi*u(1))+0.3*sin(3*pi*u(1))+0.1*sin(5*pi*u(1));
for k=2:1002
    u(k)=sin(2*pi*k/200);
    g(k)=0.6*sin(pi*u(k))+0.3*sin(3*pi*u(k))+0.1*sin(5*pi*u(k));
    y(k+1)=0.3*y(k)+0.6*y(k-1)+g(k);
end
plot(y)
legend('实际图像')

M=10;   %设置模糊系统规则个数（模糊规则是3输入1输出的）
I=3;
J=1;
% xtrain=[y(2:201) y(1:200) g(2:201)'];   %训练模糊系统参数的数据
% ytrain=y(3:202);

xtrain=[y(2:601) y(1:600) g(2:601)'];   %训练模糊系统参数的数据
ytrain=y(3:602);
%% 下面几部分实现Frank-Wolfe算法
tic%计算F-W算法所用的时间
epsilon=10e-6;
x0=10*rand(1,M*(2*I+J));%设置初始点 
x0_save=x0;%保存初始点
Tmax=100;%设置最大迭代次数

rho=0.5;%%%%%%%
alpha=0.4;%%%%% 两个是Armijo非线性线搜索的参数


%% 设置求解线性规划时变量的上下界
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
while t<Tmax %限制迭代次数为**
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
    while m< 20 %非线性规划的Armijo方法
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
legend('损失(目标)函数')

x_predict=[y(2:1001) y(1:1000) g(2:1001)'];   %训练模糊系统参数的数据
y_predict=y(3:1002);
[ df, Err,fls_f,mu_df,c_df,sigma_df] = compute_df_f( x_predict,...
                      y_predict,x0,M,I,J );%此处调用compute_df_f函数的原因：
figure                                     %    compute_df_f函数中有部分是计算模糊系统的输出fls_f，需要使用fls_f
plot(fls_f)
hold on
plot(y)
legend('预测','实际图像')
title('用F-W算法优化模糊系统参数后的预测输出与实际图像')
figure
plot(fls_f)
legend('预测')
