close all
clear
clc

%% 生成数据
% tao=31;
% Dt=tao;
% N=40020;%设置点数N，Mackey_Glass默认时隙间隔为0.1;
% [y,t]=Mackey_Glass(N,tao);%龙格-库塔(Runge-Kutta)方法求解Mackey Glass方程
% 
% n_train = 4000;
% x_star = zeros(1,n_train);
% for i = 1:n_train
%     x_star(i) = y(1+i*10);
% end
%% 
epsilon=10e-6;          %设置FW算法跳出循环时的梯度范数误差
error_precision=0.05;      %训练精度
alpha_BP=0.1;%0.05MG方程
Tmax=1000;
gama=0.00;eta=0.0005;


MC_num=1;
I = 2;%模糊规则前件个数为4
n = 2;
J = 1;%模糊规则后件数为1
data_size1 = 504;  %生成的训练数据是1000(1004-4)组
data_size2 = 100;
datanum=1;
%% 只用一定规模的数据训练模糊逻辑系统，并预测后500组的输出，例如：用1：500训练，
% 预测501：1000，...,用1：3500训练，预测3501：4000
% xtrain = [x_star(1:datanum*data_size1-4);x_star(2:datanum*data_size1-3);x_star(3:datanum*data_size1-2);x_star(4:datanum*data_size1-1);]';
% ytrain = x_star(5:datanum*data_size1)';
% x_predict = [x_star(datanum*data_size2-3:(datanum+1)*data_size2-4);x_star(datanum*data_size2-2:(datanum+1)*data_size2-3);...
%          x_star(datanum*data_size2-1:(datanum+1)*data_size2-2);x_star(datanum*data_size2:(datanum+1)*data_size2-1)]';
% y_predict = x_star(datanum*data_size2+1:(datanum+1)*data_size2)';


%%
y(1) = 0;
y(2) = 0.05;
G(1) = 0;
r = @(x)sin(2*pi*x/25)+sin(pi*x/50)+sin(2*pi*x/250); 
g = @(x,y)(x*y*(x+0.25))/(1+x^2+y^2);
K = 601;  %数据对数
for i1 = 2:K
    y(i1+1) = 0.6*y(i1)+0.2*y(i1-1)+r(i1);
end
for j = 2:K
    G(j) = g(y(j),y(j-1));
end
xtrain = [y(2:500)' y(1:500-1)'];
ytrain = G(2:500)';
x_predict = [y(502:601)' y(501:600)'];
y_predict = G(502:601)';



RMSE_SAVE=zeros(4,Tmax,10,MC_num);

for MC=1:MC_num

for k=1:8
    M=5*k
    %% 利用FW算法训练参数
    % 设置求解线性规划时变量的上下界
    c1=0*ones(1,M);
    c2=1.5*ones(1,M);
    mu1=0*ones(1,I*M);
    mu2=2*ones(1,I*M);
    sigma_ini1=0.001*ones(1,I*M);%注意sigma初始化时不可为0向量
    sigma_ini2=4*ones(1,I*M);
    bound_low=[c1,mu1, sigma_ini1]';
    bound_up=[c2,mu2, sigma_ini2]';
    for i = 1:length(bound_up)
        x0(i) = bound_up(i)*rand(1);
    end 

%     [xt_opt_FW,Error_train_FW,Error_FW,time_train_FW,t_FW] = Train_T1M_FW(xtrain,ytrain,x0,bound_low,bound_up,M,I,error_precision,epsilon,Tmax);
%     RMSE_SAVE(1,:,k,MC)=Error_FW;
%     [xt_opt_awaystepFW,Error_AFW,Error_train_awaystepFW,time_train_awaystepFW,t_awaystepFW] = Train_T1M_awaystepFW(xtrain,ytrain,x0,bound_low,bound_up,M,I,error_precision,epsilon,Tmax);
%     RMSE_SAVE(3,:,k,MC)=Error_train_awaystepFW;
    %% 利用BP算法训练参数
    % 参数初始化
%     meanF_BP=1*rand(M,I);
%     stdF_BP=1*ones(M,I);
%     yl_BP=1.5*rand(M,1);

    meanF_BP=24*rand(M,I)+12;
    stdF_BP=20*rand(M,I);
    yl_BP=10*rand(M,1);


    err_BP1 = 2000;
    t_BP=0;
    tic
    while err_BP1 > error_precision && t_BP<=Tmax
        %[meanF_BP,stdF_BP,yl_BP]=Train_T1M_BP(xtrain,ytrain,meanF,stdF,yl,alpha_BP)
        [meanF_BP,stdF_BP,yl_BP]=train_sfls_type1(xtrain,ytrain,meanF_BP,stdF_BP,yl_BP,alpha_BP);
        fls_BP1 = ST1M(xtrain,meanF_BP,stdF_BP,yl_BP);
        err_BP1=sqrt((sum(fls_BP1-ytrain).^2)/length(fls_BP1))
        RMSE_SAVE(1,t_BP+1,k,MC)=err_BP1;
        t_BP=t_BP+1
        
    end
    BP_time=toc
    %% 利用动量法
%     meanF_momentum1=0.5*rand(M,I);
%     stdF_momentum1=0.5*ones(M,I);
%     yl_momentum1=1.5*rand(M,1);
%     err_momentum = 20;
%     t_momentum=0;
%     meanF1_v=zeros(M,n);
%     stdF1_v=zeros(M,n);
%     yl_v=zeros(M,J);
%     
% 
%     tic
%     while err_momentum > error_precision && t_momentum<Tmax
%         %stdF1_v
%         meanF1_v;
%         [meanF_momentum,stdF_momentum,yl_momentum,stdF1_v,meanF1_v,yl_v]=Train_T1M_momentum(xtrain,...
%             ytrain,meanF_momentum1,stdF_momentum1,yl_momentum1,gama,eta,stdF1_v,meanF1_v,yl_v);
%         fls_momentum1 = ST1M(xtrain,meanF_momentum,stdF_momentum,yl_momentum);
%         err_momentum1=sqrt((sum(fls_momentum1-ytrain).^2)/length(fls_momentum1))
%         RMSE_SAVE(4,t_momentum+1,k,MC)=err_momentum1;
%  
%         t_momentum=t_momentum+1
%     end
%     momentum_time=toc   
    
    %% 预测
    % FW算法 
%     x_FW=xt_opt_FW;
%     c=x_FW(1:M)'; %对模糊系统的参数进行赋值
%     mu=reshape(x_FW(M+1:M+M*I),I,[]);
%     mu=mu';%
%     sigma=reshape(x_FW(M*I+M+1:M*I*2+M),I,[]);
%     sigma=sigma';
%     fls_FW = ST1M(x_predict,mu,sigma,c);
%     err_FW=sqrt((sum(fls_FW-y_predict).^2)/length(fls_FW));

%     %awaystepFW算法
%     x_awaystepFW=xt_opt_awaystepFW;
%     c=x_awaystepFW(1:M)'; %对模糊系统的参数进行赋值
%     mu=reshape(x_awaystepFW(M+1:M+M*I),I,[]);
%     mu=mu';%
%     sigma=reshape(x_awaystepFW(M*I+M+1:M*I*2+M),I,[]);
%     sigma=sigma';
%     fls_awaystepFW = ST1M(x_predict,mu,sigma,c);
%     err_awaystepFW=sqrt((sum(fls_awaystepFW-y_predict).^2)/length(fls_awaystepFW));


    % BP算法 
    fls_BP = ST1M(x_predict,meanF_BP,stdF_BP,yl_BP);
    err_BP=sqrt((sum(fls_BP-y_predict).^2)/length(fls_BP));

    % 动量
%     fls_momentum = ST1M(x_predict,meanF_momentum,stdF_momentum,yl_momentum);
%     err_momentum=sqrt((sum(fls_momentum-y_predict).^2)/length(fls_momentum));

%     table(k,:,MC)=[t_FW,Error_train_FW,err_FW,time_train_FW,...
%         t_BP,err_BP1,err_BP,BP_time];
%     table(k,:,MC)=[t_BP,err_BP1,err_BP,BP_time,...
%         t_momentum,err_momentum1,err_momentum,momentum_time]; 
      table(k,:,MC)=[t_BP,err_BP1,err_BP,BP_time];   
%      table(k,:,MC)=[t_awaystepFW,Error_AFW,err_awaystepFW,time_train_awaystepFW]; 
%       table(k,:,MC)=[t_momentum,err_momentum1,err_momentum,momentum_time];
    
    
    k=k+1;
end  %规则数量的循环

end  %蒙特卡罗次数


Table_record=table/MC_num;
save('AFW_FW_BP_dataSAVE','table','Table_record','RMSE_SAVE')


figure
plot(RMSE_SAVE(1,:,1))
hold on 
plot(RMSE_SAVE(1,:,2))
plot(RMSE_SAVE(1,:,3))
plot(RMSE_SAVE(1,:,4))
plot(RMSE_SAVE(1,:,5))
plot(RMSE_SAVE(1,:,6))
plot(RMSE_SAVE(1,:,7))
plot(RMSE_SAVE(1,:,8))
% hold off 
% figure
% plot(RMSE_SAVE(2,:,1))
% figure
% plot(RMSE_SAVE(3,:,1))
% hold on 
% plot(RMSE_SAVE(3,:,2))
% plot(RMSE_SAVE(3,:,3))
% plot(RMSE_SAVE(3,:,4))
% plot(RMSE_SAVE(3,:,5))
% plot(RMSE_SAVE(3,:,6))
% plot(RMSE_SAVE(3,:,7))
% plot(RMSE_SAVE(3,:,8))
% hold off

% figure
% plot(RMSE_SAVE(4,:,1))
% hold on 
% plot(RMSE_SAVE(4,:,2))
% plot(RMSE_SAVE(4,:,3))
% plot(RMSE_SAVE(4,:,4))
% plot(RMSE_SAVE(4,:,5))
% plot(RMSE_SAVE(4,:,6))
% plot(RMSE_SAVE(4,:,7))
% plot(RMSE_SAVE(4,:,8))
% hold off









