%用Frank-Wolfe算法调节模糊系统，使其辨识差分方程y(k+1)=0.3*y(k)+0.6*y(k-1)+g[u(k)],
%其中g(u)=0.6sin(pi*u)+0.3sin(3*pi*u)+0.1sin(5*pi*u),u(k)=cos(2*pi*k/25)
%模糊系统是3输入1输出
%% 
% clear
% close all
% clc
% 
% %% 产生数据
% y=rand(2,1);
% u(1)=sin(2*pi/200);
% g(1)=0.6*sin(pi*u(1))+0.3*sin(3*pi*u(1))+0.1*sin(5*pi*u(1));
% for k=2:3100
%     u(k)=sin(2*pi*k/200);
%     g(k)=0.6*sin(pi*u(k))+0.3*sin(3*pi*u(k))+0.1*sin(5*pi*u(k));
%     y(k+1)=0.3*y(k)+0.6*y(k-1)+g(k);
% end
% 
% % I=3;    %前件个数
% I=3;    %前件个数
% J=1;    %后件个数
% con=0; %数据起始位置
% 
% train_num=600;
% test_num=400;
% xtrain=[y(con+I-1:con+train_num+I-2) y(con+I-2:con+train_num+I-3) g(con+I-1:con+train_num+I-2)']; %训练模糊系统参数的数据
% ytrain=y(con+I:con+train_num+I-1);
% 
% x_predict=[y(con+train_num+I-1:con+train_num+test_num+I-2) y(con+train_num+I-2:...
%     con+train_num+test_num+I-3) g(con+train_num+I-1:con+train_num+test_num+I-2)'];   %预测模糊系统参数的数据
% y_predict=y(con+train_num+I:con+train_num+test_num+I-1);
% 
number=50;
% table=zeros(number/5,9);
%% 
% 参数初始化
epsilon=10e-6;          %设置FW算法跳出循环时的梯度范数误差
error_precision =1;      %训练精度
nnn=5:5:number;
%% 
for z=1:length(nnn)
z
M=nnn(z);   %设置模糊系统规则个数（模糊规则是3输入1输出的）可以修改
%%
% t1=100;
% xtrain=rand(t1,I);
% ytrain=2*rand(t1,1);
% x_predict=rand(t1/2,I);
% y_predict=rand(t1/2,1);
%% 利用FW算法训练参数
% 设置求解线性规划时变量的上下界
c_mu1=-5*ones(1,I*M);
c_mu2=5*ones(1,I*M);
for i=1:I*M
    if mod(i,3)==0
       c_mu1(i+M)=-1; 
       c_mu2(i+M)=1;
    end
end
sigma_ini1=0.01*ones(1,I*M);%注意sigma初始化时不可为0向量
sigma_ini2=4*ones(1,I*M);
bound_low=[c_mu1 sigma_ini1]';
bound_up=[c_mu2 sigma_ini2]';
for i = 1:length(bound_up)
    x0(i) = bound_up(i)*rand(1);
end 

[xt_opt_FW,Error_train_FW,time_train_FW,t_FW] = Train_T1M_FW(xtrain,ytrain,x0,bound_low,bound_up,M,I,error_precision,epsilon);

[xt_opt_awaystepFW,Error_train_awaystepFW,time_train_awaystepFW,t_awaystepFW] = Train_T1M_awaystepFW(xtrain,ytrain,x0,bound_low,bound_up,M,I,error_precision,epsilon);
%% 利用BP算法训练参数
% 参数初始化
meanF=-10+5*rand(M,I);
stdF=4*ones(M,I);
yl=-10+5*rand(M,1);
alpha=2;
%epoch=length(alpha);
err_BP1 = 20;
t_BP=0;
tic
V=0;S=0;
while err_BP1 > error_precision && t_BP<=5000
    [meanF,stdF,yl,V,S]=Train_T1M_BP(xtrain,ytrain,meanF,stdF,yl,alpha,V,S);
    fls_BP1 = ST1M(xtrain,meanF,stdF,yl);
    err_BP1=0.5*(norm(fls_BP1-ytrain)^2);
    t_BP=t_BP+1
end
BP_time=toc
%% 预测
% FW算法 
x_FW=xt_opt_FW;
c=x_FW(1:M)'; %对模糊系统的参数进行赋值
mu=reshape(x_FW(M+1:M+M*I),I,[]);
mu=mu';%
sigma=reshape(x_FW(M*I+M+1:M*I*2+M),I,[]);
sigma=sigma';
fls_FW = ST1M(x_predict,mu,sigma,c);
err_FW=0.5*(norm(fls_FW-y_predict)^2)

%awaystepFW算法
x_awaystepFW=xt_opt_awaystepFW;
c=x_awaystepFW(1:M)'; %对模糊系统的参数进行赋值
mu=reshape(x_awaystepFW(M+1:M+M*I),I,[]);
mu=mu';%
sigma=reshape(x_awaystepFW(M*I+M+1:M*I*2+M),I,[]);
sigma=sigma';
fls_awaystepFW = ST1M(x_predict,mu,sigma,c);
err_awaystepFW=0.5*(norm(fls_awaystepFW-y_predict)^2)


% BP算法 
fls_BP = ST1M(x_predict,meanF,stdF,yl);
err_BP=0.5*(norm(fls_BP-y_predict)^2)

table(z,:)=[t_FW,err_FW,time_train_FW,t_awaystepFW,err_awaystepFW,time_train_awaystepFW,t_BP,err_BP,BP_time];
z=z+1
end
%% 画图

% figure(1)
% t=con+train_num+I:con+train_num+test_num+I-1;
% plot(t,y_predict,'k',t,fls_FW,'r',t,fls_BP,'b')
% legend('实际值','FW','BP')
% xlabel('t');
% ylabel('y(t)');
% title('用FW,BP算法优化模糊系统参数后的预测输出与实际输出对比')

