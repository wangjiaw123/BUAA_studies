
close all
clear
clc

%% 生成数据
%star_parallel_time=clock;
tao=31;
Dt=tao;
N=40020;%设置点数N，Mackey_Glass默认时隙间隔为0.1;
[y,t]=Mackey_Glass(N,tao);%龙格-库塔(Runge-Kutta)方法求解Mackey Glass方程
%% 原始图像
% figure(1)
% subplot(1,2,1)% 序列直出 
% plot(t,y,'LineWidth',1.0);
% xlabel('t')
% ylabel('s(t)')
% title('序列直出 ');
% 
% subplot(1,2,2)% 相图时差 
% plot(y(Dt*10+10000+1:end),y(10001:end-Dt*10),'LineWidth',0.5)
% xlabel('s(t)')
% ylabel('s(t-\tau)')
% title('相图时差');

n_train = 4000;
x_star = zeros(1,n_train);
for i = 1:n_train
    x_star(i) = y(1+i*10);
end


%N = 40;%设置规则数目20
%% 
max_select_train = 3500;
data_size = 500;

%% 全局变量
n = 4;%模糊规则前件个数为4
alpha =0.4;%KM时取0.4
setRMSE_accept=0.04;  %设置可接受的误差
%% 此处是最终要保存的数值矩阵，多维数组一开始大小不要分配错了！
RULE_NUMBER = 4;%10
RULE_WIDTH = 4;%10
CORE_NUMBER =24;%36      %查看并行部分，核心数的使用是递减的！
%CORE_NUMBER_WIDTH=1;
DATA_SIZE = 7;%7  %max_select_train/data_size=7                              
Tmax = 30;  %160                                       %设置模型训练中的最大迭代次数

Zero1 = zeros(CORE_NUMBER,DATA_SIZE,RULE_NUMBER);  %1维是分配的核心数的个数，2维是数据规模数的个数，3维是规则数的个数
SAVEparallel_total_time_A    = Zero1;              %分别是并行总时间，并行内部时间，串行时间，计算加速比、并行效率、串行比例数组的初始化
SAVEparallel_interior_time_B = Zero1;         
SAVEspeedUP_ratio_D          = Zero1;
SAVEparallel_efficiency_E    = Zero1;
SAVEserial_proportion_H     = Zero1;

SAVE_RMSE_parallel_G = zeros(CORE_NUMBER,Tmax,DATA_SIZE,RULE_NUMBER);   %1维是分配的核心数的个数，2维是Tmax(这一维用于存储RMSE),3维是数据规模数的个数，4维是规则数的个数

SAVE_parallel_M1    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_M2    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_sigma = zeros(CORE_NUMBER,RULE_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_c1    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_WIDTH,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_c2    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_WIDTH,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_Ytest = zeros(CORE_NUMBER,max_select_train,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_Ypridect = zeros(CORE_NUMBER,data_size,DATA_SIZE,RULE_NUMBER);

SAVEserial_total_time_C = zeros(RULE_NUMBER,DATA_SIZE);        %1维是规则数的个数，2维是数据规模数的个数
SAVE_RMSE_serial_F    = zeros(RULE_NUMBER,Tmax,DATA_SIZE);     %维是规则数的个数，2维是Tmax(这一维用于存储RMSE),3维是数据规模数的个数

SAVE_serial_M1    = zeros(RULE_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULE_NUMBER);
SAVE_serial_M2    = zeros(RULE_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULE_NUMBER);
SAVE_serial_sigma = zeros(RULE_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULE_NUMBER);

SAVE_serial_c1    = zeros(RULE_NUMBER*RULE_WIDTH,DATA_SIZE,RULE_NUMBER);
SAVE_serial_c2    = zeros(RULE_NUMBER*RULE_WIDTH,DATA_SIZE,RULE_NUMBER);
SAVE_serial_Ytest = zeros(max_select_train,DATA_SIZE,RULE_NUMBER);
SAVE_serial_Ypridect = zeros(data_size,DATA_SIZE,RULE_NUMBER);         



%% 串行程序部分

for rule_num=2:RULE_NUMBER
    
    N=RULE_WIDTH^rule_num;          %规则数目的循环
    
for i=1:DATA_SIZE                   %任务规模的循环
   %% 时间记录参数初始化
    serial_total_time = 0;   
    %% 选择训练数据与测试数据
    % 只用一定规模的数据训练模糊逻辑系统，并预测后500组的输出，例如：用1：500训练，
    % 预测501：1000，...,用1：3500训练，预测3501：4000
    X_train = [x_star(1:i*data_size-4);x_star(2:i*data_size-3);x_star(3:i*data_size-2);x_star(4:i*data_size-1);]';
    Y_train = x_star(5:i*data_size)';
    X_test = [x_star(i*data_size-3:(i+1)*data_size-4);x_star(i*data_size-2:(i+1)*data_size-3);...
            x_star(i*data_size-1:(i+1)*data_size-2);x_star(i*data_size:(i+1)*data_size-1)]';
    Y_test = x_star(i*data_size+1:(i+1)*data_size)';
    
    %% 训练参数初始化
    M1 = rand(N,n);
    M2 = M1 + 0.8;
    sigma = rand(N,n);
    c1 = rand(N,1);
    c2 = c1+1;    
    
    %% 串行方法
    fprintf('规则数为：%d,开始第 %d 种规模的训练(采用串行程序)，训练数据规模为:%d ,预测数据规模为:%d \n',N,i,data_size*i,data_size);
    star_serial_time=clock;
    R=[];
    for i1=1:length(X_train(:,1))
           [R1,R2,R(i1)]=sfls_type2(X_train(i1,:),M1,M2,sigma,c1,c2); %R(i)=1/2(R1+R2)
    end
    RMSE = sqrt(sum((R'-Y_train).^2)/length(R))
    RMSE_serial_save=zeros(1,Tmax);
    T1=0;
    M1_serial_save=zeros(N,n,Tmax);
    M2_serial_save=zeros(N,n,Tmax);
    sigma_serial_save=zeros(N,n,Tmax);
    c1_serial_save=zeros(N,1,Tmax);
    c2_serial_save=zeros(N,1,Tmax);    
    
    while  T1 <Tmax %&& RMSE  >=setRMSE_accept
           [M1,M2,c1,c2,sigma,I2l,I2u,I1u,I1l]=train_sfls_type2(X_train,...
               Y_train,M1,M2,sigma,c1,c2,alpha);
           
           for i2=1:length(X_train(:,1))
                  [R1,R2,R(i2)]=sfls_type2(X_train(i2,:),M1,M2,sigma,c1,c2);
           end   
           
           RMSE = sqrt(sum((R'-Y_train).^2)/length(R))
           if isnan(RMSE) == 1
               RMSE = 10*setRMSE_accept
               M1 = rand(N,n);
               M2 = M1 + 0.8;
               sigma = rand(N,n);
               c1 = rand(N,1);
               c2 = c1+1; 
           end
           RMSE_serial_save(T1+1)=RMSE; 
           
           M1_serial_save(:,:,T1+1)=M1;
           M2_serial_save(:,:,T1+1)=M2;
           sigma_serial_save(:,:,T1+1)=sigma;
           c1_serial_save(:,:,T1+1)=c1;
           c2_serial_save(:,:,T1+1)=c2;    
           
           T1=T1+1
           
%            figure(2)
%            plot(1:length(R),[R;Y_train'],'LineWidth',1.0)
%            title(['规则规模为',num2str(rule_num),',任务规模数为',num2str(i),...
%                '，串行方法第',num2str(T1),'次训练实际值与预测值图像','(\tau =',num2str(tao),')'])
%            %xticks([data_size*(i-1)+1:fix(length(R)/5):data_size*(i-1)+length(R)]) 
%            xlabel('t')
%            ylabel('s(t)')
%            legend('预测值','实际值')      
                        
    end    

    wz_serial=find(RMSE_serial_save==min(RMSE_serial_save(RMSE_serial_save~=0)));
    wz_serial=wz_serial(1);
    Rs_serial=zeros(1,length(X_test(:,1)));   
    for ik=1:length(X_test(:,1))
            [R11,R22,Rs_serial(ik)]=sfls_type2(X_test(ik,:),M1_serial_save(:,:,wz_serial),...
                M2_serial_save(:,:,wz_serial),sigma_serial_save(:,:,wz_serial),...
                c1_serial_save(:,:,wz_serial),c2_serial_save(:,:,wz_serial));
    end
    RMSE = sqrt(sum((Rs_serial'-Y_test).^2)/length(Rs_serial))
    
    end_serial_time=clock;
    serial_total_time=etime(end_serial_time,star_serial_time)  
    SAVEserial_total_time_C(rule_num,i)      = serial_total_time;   
    
    SAVE_serial_M1(1:length(M1_serial_save(:,1,wz_serial)),:,i,rule_num)    = M1_serial_save(:,:,wz_serial);
    SAVE_serial_M2(1:length(M2_serial_save(:,1,wz_serial)),:,i,rule_num)    = M2_serial_save(:,:,wz_serial);
    SAVE_serial_sigma(1:length(sigma_serial_save(:,1,wz_serial)),:,i,rule_num) = sigma_serial_save(:,:,wz_serial);
    SAVE_serial_c1(1:length(c1_serial_save(:,1,wz_serial)),i,rule_num)      = c1_serial_save(:,:,wz_serial);
    SAVE_serial_c2(1:length(c2_serial_save(:,1,wz_serial)),i,rule_num)      = c2_serial_save(:,:,wz_serial);
    SAVE_serial_Ytest(1:length(Y_test),i,rule_num)       = Y_test;
    SAVE_serial_Ypridect(1:length(Rs_serial),i,rule_num) = Rs_serial;
    SAVE_RMSE_serial_F(rule_num,:,i)    = RMSE_serial_save; 
%     figure(3)
%     subplot(1,2,1)
%     plot(i*data_size+1:i*data_size+length(Rs_serial),Rs_serial,'LineWidth',1.0)
%     xlabel('t')
%     ylabel('s(t)')
%     title('序列直出 ');
%     subplot(1,2,2)
%     plot(Rs_serial(tao+1:end),Rs_serial(1:end-tao),'LineWidth',0.5)
%     xlabel('s(t)')
%     ylabel('s(t-\tau)')
%     title('相图时差');

%     figure(4)
%     plot(i*data_size+1+length(R):i*data_size+length(Rs_serial)+length(R),[Rs_serial;Y_test'],'LineWidth',1.0)
%     xlabel('t')
%     ylabel('s(t)')
%     title(['串行方法实际值与最终预测值图像','(\tau =',num2str(tao),')'])
%     legend('预测值','实际值')    
  
end
        
end

%% 并行程序部分
for core_num=2:2:CORE_NUMBER %   %核心数的循环
    
    parpool(core_num);          %在MATLAB中开启核心数为core_num的并行计算池  
    
for rule_num=2:RULE_NUMBER
    
    N=RULE_WIDTH^rule_num;                     %规则数目的循环
    %N=RULE(rule_num); 
for i=1:DATA_SIZE                      %任务规模的循环
   %% 时间记录参数初始化
    parallel_total_time = 0;
    parallel_interior_time = 0;
    %% 选择训练数据与测试数据
    % 只用一定规模的数据训练模糊逻辑系统，并预测后500组的输出，例如：用1：500训练，
    % 预测501：1000，...,用3001：3500训练，预测3501：4000    X_train = [x_star(1:i*data_size-4);x_star(2:i*data_size-3);x_star(3:i*data_size-2);x_star(4:i*data_size-1);]';
    X_train = [x_star(1:i*data_size-4);x_star(2:i*data_size-3);x_star(3:i*data_size-2);x_star(4:i*data_size-1);]';
    Y_train = x_star(5:i*data_size)';
    X_test = [x_star(i*data_size-3:(i+1)*data_size-4);x_star(i*data_size-2:(i+1)*data_size-3);...
            x_star(i*data_size-1:(i+1)*data_size-2);x_star(i*data_size:(i+1)*data_size-1)]';
    Y_test = x_star(i*data_size+1:(i+1)*data_size)';
     
    %% 训练参数初始化
    M1 = rand(N,n);
    M2 = M1 + 1;%0.3
    sigma = rand(N,n);
    c1 = rand(N,1);
    c2 = c1+1.5; %0.55
    %% 并行方法
    star_parallel_time=clock;%直到程序结束，计算并行程序总的运行时间
    fprintf('规则数为：%d ,开始第 %d 种规模的训练(采用并行程序)，训练数据规模为:%d ,预测数据规模为:%d \n',N,i,data_size*i,data_size)
    R=[];
    % 第1次计算2型模糊系统的输出
    tic
    parfor i11=1:length(X_train(:,1))
           [R1,R2,R(i11)]=sfls_type2(X_train(i11,:),M1,M2,sigma,c1,c2); %R(i)=1/2(R1+R2)
    end
    time11=toc;
    parallel_interior_time=parallel_interior_time+time11;

    RMSE = sqrt(sum((R'-Y_train).^2)/length(R))
    % 使用并行算法训练模型
    % 设置参数的储存变量
    M1_save=zeros(N,n,Tmax);
    M2_save=zeros(N,n,Tmax);
    c1_save=zeros(N,1,Tmax);
    c2_save=zeros(N,1,Tmax);
    sigma_save=zeros(N,n,Tmax);
    RMSE_parallel_save=zeros(1,Tmax);
    
    row_num=length(X_train(:,1));%获取当前训练数据对的组数
    lack_num=core_num-mod(row_num,core_num);%计算填充数
    R_save=zeros(Tmax,length(Y_train)+lack_num);
    for k=1:lack_num
        X_train=[X_train' X_train(k,:)']';
        Y_train=[Y_train' Y_train(k,:)']';
    end
    total_num=row_num+lack_num;
    circulation_num=total_num/core_num; %分配每个核心上的训练数据规模
    XX_train=zeros(core_num,length(X_train));
    
    time21=0;
    time22=0;
    T2=0;
    %% 使用Adam优化方法之前的参数设置，内存分配
    beta1 = 0.9;beta2=0.999;epsilon=10e-6;
    m_g_sigma = zeros(size(sigma));
    m_g_M1 = zeros(size(M1));
    m_g_M2 = zeros(size(M2));
    m_g_c1 = zeros(size(c1));
    m_g_c2 = zeros(size(c2));
    
    v_g_sigma = zeros(size(sigma));
    v_g_M1 = zeros(size(M1));
    v_g_M2 = zeros(size(M2));
    v_g_c1 = zeros(size(c1));
    v_g_c2 = zeros(size(c2)); 
    ONoff_parallel_flag=0;
    %%     
    while  T2 <Tmax %&& RMSE >=setRMSE_accept
        
           sigma_help=zeros(N,n,core_num);
           M1_help=zeros(N,n,core_num);
           M2_help=zeros(N,n,core_num);
           c1_help=zeros(N,1,core_num);
           c2_help=zeros(N,1,core_num);  
           sigma22=sigma;M22=M2;M12=M1;
           c12=c1;c22=c2;   
           XX_train=zeros(circulation_num,n,core_num);
           YY_train=zeros(circulation_num,1,core_num); 
           
          %% 训练数据随机重排
           X_copy=zeros(size(X_train));
           Y_copy=zeros(size(Y_train));
           lengthX=length(X_train(:,1));
           RandID=randperm(lengthX);     
           for ji = 1:lengthX
               X_copy(ji,:) = X_train(RandID(ji),:);
               Y_copy(ji,:) = Y_train(RandID(ji),:);
           end
           %% 
           tic
           parfor i22=1:core_num
                  XX_train(:,:,i22)=X_copy(((i22-1)*circulation_num+1):(i22*circulation_num),:);
                  YY_train(:,:,i22)=Y_copy(((i22-1)*circulation_num+1):(i22*circulation_num));
                  [M1_help(:,:,i22),M2_help(:,:,i22),c1_help(:,:,i22),c2_help(:,:,i22),sigma_help(:,:,i22)]...
                            =train_sfls_type2(XX_train(:,:,i22),YY_train(:,:,i22),M1,M2,sigma,c1,c2,alpha*core_num);%2
           end
           ttime21=toc;
           time21=time21+ttime21;
           
           for ii=2:core_num
               M1_help(:,:,1)=M1_help(:,:,1)+M1_help(:,:,ii);
               M2_help(:,:,1)=M2_help(:,:,1)+M2_help(:,:,ii);
               sigma_help(:,:,1)=sigma_help(:,:,1)+sigma_help(:,:,ii);
               c1_help(:,:,1)=c1_help(:,:,1)+c1_help(:,:,ii);
               c2_help(:,:,1)=c2_help(:,:,1)+c2_help(:,:,ii);
           end     

%            sigma=sigma_help(:,:,1)/core_num;
%            M1=M1_help(:,:,1)/core_num;
%            M2=M2_help(:,:,1)/core_num;
%            c1=c1_help(:,:,1)/core_num;
%            c2=c2_help(:,:,1)/core_num;   
           
           g_sigma=sigma_help(:,:,1)/core_num;
           g_M1=M1_help(:,:,1)/core_num;
           g_M2=M2_help(:,:,1)/core_num;
           g_c1=c1_help(:,:,1)/core_num;
           g_c2=c2_help(:,:,1)/core_num;  
           
%            sigma = g_sigma;
%            M1 = g_M1;
%            M2 = g_M2;
%            c1 = g_c1;
%            c2 = g_c2;           
          %% 实现Adam算法
           m_g_sigma = beta1 * m_g_sigma + (1-beta1)* g_sigma;
           m_g_M1 = beta1 * m_g_M1 + (1-beta1) * g_M1;
           m_g_M2 = beta1 * m_g_M2 + (1-beta1) * g_M2;
           m_g_c1 = beta1 * m_g_c1 + (1-beta1) * g_c1;
           m_g_c2 = beta1 * m_g_c2 + (1-beta1) * g_c2;
   
           v_g_sigma = beta2 * v_g_sigma + (1-beta2) * (g_sigma.*g_sigma);
           v_g_M1 = beta2 * v_g_M1 + (1-beta2) * (g_M1.*g_M1);
           v_g_M2 = beta2 * v_g_M2 + (1-beta2) * (g_M2.*g_M2);
           v_g_c1 = beta2 * v_g_c1 + (1-beta2) * (g_c1.*g_c1);
           v_g_c2 = beta2 * v_g_c2 + (1-beta2) * (g_c2.*g_c2);              
           
           m_g_sigma = m_g_sigma/(1-beta1^(T2+1));
           m_g_M1 = m_g_M1/(1-beta1^(T2+1));
           m_g_M2 = m_g_M2/(1-beta1^(T2+1));
           m_g_c1 = m_g_c1/(1-beta1^(T2+1));
           m_g_c2 = m_g_c2/(1-beta1^(T2+1));
           
           v_g_sigma = v_g_sigma/(1-beta2^(T2+1));
           v_g_M1 = v_g_M1/(1-beta2^(T2+1));
           v_g_M2 = v_g_M2/(1-beta2^(T2+1));
           v_g_c1 = v_g_c1/(1-beta2^(T2+1));
           v_g_c2 = v_g_c2/(1-beta2^(T2+1));          
           alpha1=0.02;
           sigma = g_sigma - alpha1 * (m_g_sigma./(sqrt(v_g_sigma)+epsilon));
           M1 = g_M1 - alpha1 * (m_g_M1./(sqrt(v_g_M1)+epsilon));
           M2 = g_M2 - alpha1 * (m_g_M2./(sqrt(v_g_M2)+epsilon));
           c1 = g_c1 - alpha1 * (m_g_c1./(sqrt(v_g_c1)+epsilon));
           c2 = g_c2 - alpha1 * (m_g_c2./(sqrt(v_g_c2)+epsilon));
           
          %%           
           M1_save(:,:,T2+1)=M1;
           M2_save(:,:,T2+1)=M2;
           c1_save(:,:,T2+1)=c1;
           c2_save(:,:,T2+1)=c2;
           sigma_save(:,:,T2+1)=sigma;
      
           tic
           parfor i33=1:length(X_train(:,1))
                  [R1,R2,R(i33)]=sfls_type2(X_train(i33,:),M1,M2,sigma,c1,c2);
           end
           ttime22=toc;
           time22=time22+ttime22;
     
           RMSE = sqrt(sum((R'-Y_train).^2)/length(R))
           if isnan(RMSE) == 1
               RMSE = 10*setRMSE_accept
               M1 = rand(N,n);
               M2 = M1 + 1;
               sigma = rand(N,n);
               c1 = rand(N,1);
               c2 = c1+1.5; 
           end
           RMSE_parallel_save(T2+1)=RMSE;
      
           T2 = T2+1
           
%            figure(5)
%            plot(1:length(R),[R;Y_train'],'LineWidth',1.0)
%            title(['核心数为',num2str(core_num),'规则规模为',num2str(rule_num),',任务规模数为',num2str(i),...
%               '，并行方法第',num2str(T2),'次训练实际值与预测值图像','(\tau =',num2str(tao),')'])
%            %xticks([data_size*(i-1)+1:fix(length(R)/5):data_size*(i-1)+length(R)]) 
%            xlabel('t')
%            ylabel('s(t)')
%            legend('预测值','实际值')

    end
    parallel_interior_time=parallel_interior_time+time21+time22;   
    
    wz=find(RMSE_parallel_save==min(RMSE_parallel_save(RMSE_parallel_save~=0)));
    wz=wz(1);
    Rs=zeros(1,length(X_test(:,1)));
    
    tic
    parfor i=1:length(X_test(:,1))
            [R11,R22,Rs(i)]=sfls_type2(X_test(i,:),M1_save(:,:,wz),M2_save(:,:,wz),sigma_save(:,:,wz),c1_save(:,:,wz),c2_save(:,:,wz));
    end
    time12=toc;
    parallel_interior_time=parallel_interior_time+time12;
    
    RMSE = sqrt(sum((Rs'-Y_test).^2)/length(Rs))
    
    end_parallel_time=clock;
    parallel_total_time=etime(end_parallel_time,star_parallel_time)  
    
    SAVE_parallel_M1(core_num,1:length(M1_save(:,1,wz)),:,i,rule_num)    = M1_save(:,:,wz);
    SAVE_parallel_M2(core_num,1:length(M2_save(:,1,wz)),:,i,rule_num)    = M2_save(:,:,wz);
    SAVE_parallel_sigma(core_num,1:length(sigma_save(:,1,wz)),:,i,rule_num) = sigma_save(:,:,wz);
    SAVE_parallel_c1(core_num,1:length(c1_save(:,1,wz)),i,rule_num)      = c1_save(:,:,wz);
    SAVE_parallel_c2(core_num,1:length(c2_save(:,1,wz)),i,rule_num)      = c2_save(:,:,wz);
    SAVE_parallel_Ytest(core_num,1:length(Y_test),i,rule_num)       = Y_test;
    SAVE_parallel_Ypridect(core_num,1:length(Rs),i,rule_num) = Rs;   

%     figure(6)
%     subplot(1,2,1)
%     plot(i*data_size+1:i*data_size+length(Rs),Rs,'LineWidth',1.0)
%     xlabel('t')
%     ylabel('s(t)')
%     title('序列直出 ');
%     subplot(1,2,2)
%     plot(Rs(tao+1:end),Rs(1:end-tao),'LineWidth',0.5)
%     xlabel('s(t)')
%     ylabel('s(t-\tau)')
%     title('相图时差');
% 
%     figure(7)
%     plot(i*data_size+1+length(R):i*data_size+length(Rs)+length(R),[Rs;Y_test'],'LineWidth',1.0)
%     xlabel('t')
%     ylabel('s(t)')
%     title(['并行方法实际值与最终预测值图像','(\tau =',num2str(tao),')'])
%     legend('预测值','实际值')

    %% 记录并行总时间，并行内部时间，串行时间，计算加速比、并行效率,串行比例
    SAVEparallel_total_time_A(core_num,i,rule_num)    = parallel_total_time;      
    SAVEparallel_interior_time_B(core_num,i,rule_num) = parallel_interior_time;         
    SAVEspeedUP_ratio_D(core_num,i,rule_num)          = SAVEserial_total_time_C(rule_num,i)/parallel_total_time;
    SAVEparallel_efficiency_E(core_num,i,rule_num)    = SAVEserial_total_time_C(rule_num,i)/(core_num*parallel_total_time);
    SAVEserial_proportion_H(core_num,i,rule_num)     = (parallel_total_time - parallel_interior_time)/parallel_total_time;   
    SAVE_RMSE_parallel_G(core_num,:,i,rule_num) = RMSE_parallel_save;
end       
    
end  
   
   delete(gcp('nocreate')) %每一轮循环后关闭MATLAB并行计算池，下一轮开启不同核心数的并行计算池
end

%% 保存所有需要的数据到SAVE_DATA.mat文件中，以便于列表或做图
save('SAVEParallel_performance_parameter.mat','SAVEparallel_total_time_A','SAVEparallel_interior_time_B',...
    'SAVEserial_total_time_C','SAVEspeedUP_ratio_D','SAVEparallel_efficiency_E',...
    'SAVE_RMSE_serial_F','SAVE_RMSE_parallel_G','SAVEserial_proportion_H')
save('SAVEdata_parallel_parameter','SAVE_parallel_M1','SAVE_parallel_M2','SAVE_parallel_sigma',...
    'SAVE_parallel_c1','SAVE_parallel_c2','SAVE_parallel_Ytest','SAVE_parallel_Ypridect')
save('SAVEdata_serial_parameter','SAVE_serial_M1','SAVE_serial_M2','SAVE_serial_sigma',...
    'SAVE_serial_c1','SAVE_serial_c2','SAVE_serial_Ytest','SAVE_serial_Ypridect')




