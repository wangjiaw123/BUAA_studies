
close all
clear
clc

%% 生成数据
y(1) = 0;
y(2) = 0.05;
G(1) = 0;
r = @(x)sin(2*pi*x/25)%+sin(pi*x/50)+sin(2*pi*x/250); 
g = @(x,y)(x*y*(x+0.25))/(1+x^2+y^2);
K = 500;  %数据对数
for i1 = 2:K
    y(i1+1) = 0.6*y(i1)+0.2*y(i1-1)+r(i1);
end
for j = 2:K
    G(j) = g(y(j),y(j-1));
end
Xtrain = [y(2:K)' y(1:K-1)'];
Ytrain = G(2:K)';


%% 参数初始化
[M,n] = size(Xtrain);
Tmax = 40;
alpha = 0.004;
setRMSE_accept = 0.04;
RULE_WIDTH = 8;
RULE_NUMBER = 8;
CORE_NUMBER = 24;

SAVEserial_total_time = zeros(RULE_NUMBER);
SAVE_serial_M1    = zeros(RULE_NUMBER*RULE_WIDTH,n,RULE_NUMBER);
SAVE_serial_M2    = zeros(RULE_NUMBER*RULE_WIDTH,n,RULE_NUMBER);
SAVE_serial_sigma = zeros(RULE_NUMBER*RULE_WIDTH,n,RULE_NUMBER);
SAVE_serial_C    = zeros(RULE_NUMBER*RULE_WIDTH,n+1,RULE_NUMBER);
SAVE_serial_S    = zeros(RULE_NUMBER*RULE_WIDTH,n+1,RULE_NUMBER);
SAVE_serial_Ytest = zeros(K,RULE_NUMBER);
SAVE_serial_Ypridect = zeros(K,RULE_NUMBER);
SAVE_RMSE_serial = zeros(Tmax,RULE_NUMBER);

DATA_SIZE=1;
Zero1 = zeros(CORE_NUMBER,DATA_SIZE,RULE_NUMBER);  %1维是分配的核心数的个数，2维是数据规模数的个数，3维是规则数的个数
SAVEparallel_total_time    = Zero1;              %分别是并行总时间，并行内部时间，串行时间，计算加速比、并行效率、串行比例数组的初始化
SAVEparallel_interior_time = Zero1;         
SAVEspeedUP_ratio          = Zero1;
SAVEparallel_efficiency    = Zero1;
SAVEserial_proportion      = Zero1;

SAVE_RMSE_parallel = zeros(CORE_NUMBER,Tmax,DATA_SIZE,RULE_NUMBER);   %1维是分配的核心数的个数，2维是Tmax(这一维用于存储RMSE),3维是数据规模数的个数，4维是规则数的个数

SAVE_parallel_M1    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_M2    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_sigma = zeros(CORE_NUMBER,RULE_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_C    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_WIDTH,n+1,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_S    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_WIDTH,n+1,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_Ytest = zeros(CORE_NUMBER,K,DATA_SIZE,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_Ypridect = zeros(CORE_NUMBER,K,DATA_SIZE,RULE_NUMBER);

%% 串行部分
for rule_num=1:RULE_NUMBER
    
    N=RULE_WIDTH*rule_num;          %规则数目的循环8,16 
    
    serial_total_time = 0;          % 时间记录参数初始化
   
    %% 训练参数初始化
    M1 = rand(N,n);
    M2 = M1 + 0.8;
    sigma = rand(N,n)+0.1;
    C = rand(N,n+1);
    S = rand(N,n+1);     
    %% 串行方法
    fprintf('采用串行程序,规则数为：%d',N);
    star_serial_time=clock;
    R=[];
    [R1,R2,R]=tsk_type2(Xtrain,M1,M2,sigma,C,S); %R(i)=1/2(R1+R2)
    RMSE = sqrt(sum((R'-Ytrain).^2)/length(R))    
    
    T1=0;
    
    RMSE_serial_save=zeros(1,Tmax);
    M1_serial_save=zeros(N,n,Tmax);
    M2_serial_save=zeros(N,n,Tmax);
    sigma_serial_save=zeros(N,n,Tmax);
    C_serial_save=zeros(N,n+1,Tmax);
    S_serial_save=zeros(N,n+1,Tmax);    
    
    while  T1 <Tmax
        
           [M1,M2,C,S,sigma]=train_tsk_type2(double(Xtrain),...
               Ytrain,M1,M2,sigma,C,S,alpha);
           
           for i2=1:length(Xtrain(:,1))
                  [R1,R2,R(i2)]=tsk_type2(Xtrain(i2,:),M1,M2,sigma,C,S);
           end    
           RMSE = sqrt(sum((R'-Ytrain).^2)/length(R))
           
           if isnan(RMSE) == 1
               RMSE = 5*setRMSE_accept
               M1 = rand(N,n);
               M2 = M1 + 0.8;
               sigma = rand(N,n);
               C = rand(N,n+1);
               S = rand(N,n+1); 
           end   
           
           RMSE_serial_save(T1+1)=RMSE;  
           M1_serial_save(:,:,T1+1)=M1;
           M2_serial_save(:,:,T1+1)=M2;
           sigma_serial_save(:,:,T1+1)=sigma;
           C_serial_save(:,:,T1+1)=C;
           S_serial_save(:,:,T1+1)=S;    
           
           T1=T1+1
           
%           figure(1)
%           plot(1:length(R),[R;Ytrain'],'LineWidth',1.0)
%           title(['规则数为',num2str(N),'，串行方法第',num2str(T1),'次训练实际值与预测值图像'])
%           xlabel('t')
%           ylabel('g(t)')
%           legend('预测值','实际值')      
                        
    end    
    wz_serial=find(RMSE_serial_save==min(RMSE_serial_save(RMSE_serial_save~=0)));
    wz_serial=wz_serial(1);
    Rs_serial=zeros(1,length(Xtrain(:,1)));   
  
    [R11,R22,Rs_serial]=tsk_type2(Xtrain,M1_serial_save(:,:,wz_serial),...
             M2_serial_save(:,:,wz_serial),sigma_serial_save(:,:,wz_serial),...
             C_serial_save(:,:,wz_serial),S_serial_save(:,:,wz_serial));
    
    RMSE = sqrt(sum((Rs_serial'-Ytrain).^2)/length(Rs_serial))
    
    end_serial_time=clock;
    serial_total_time=etime(end_serial_time,star_serial_time) 
    
    SAVEserial_total_time(rule_num)      = serial_total_time;   
    SAVE_serial_M1(1:length(M1(:,1)),:,rule_num)    = M1;
    SAVE_serial_M2(1:length(M2(:,1)),:,rule_num)     = M2;
    SAVE_serial_sigma(1:length(sigma(:,1)),:,rule_num)  = sigma;
    SAVE_serial_C(1:length(C(:,1)),:,rule_num)     = C;
    SAVE_serial_S(1:length(S(:,1)),:,rule_num)    = S;
    SAVE_serial_Yrain(1:length(Ytrain(:,1)),rule_num) = Ytrain;
    SAVE_serial_Ypridect(1:length(Rs_serial),rule_num) = Rs_serial; 
    SAVE_RMSE_serial(:,rule_num) =  RMSE_serial_save;
end    

%% 并行程序部分
for core_num=1:2:CORE_NUMBER %   %核心数的循环
    
    parpool(core_num);          %在MATLAB中开启核心数为core_num的并行计算池  
    
for rule_num=1:RULE_NUMBER
    
    N=RULE_WIDTH*rule_num;                     %规则数目的循环
    %N=RULE(rule_num); 
    DATA_SIZE=1 ;
    i=1;
   %% 时间记录参数初始化
    parallel_total_time = 0;
    parallel_interior_time = 0;
  
    %% 训练参数初始化
    M1 = rand(N,n);
    M2 = M1 + 1.5;%1
    sigma = rand(N,n)+0.2;
    C = rand(N,1+n);
    S = rand(N,1+n); %0.55
    %% 并行方法
    star_parallel_time=clock;%直到程序结束，计算并行程序总的运行时间
    fprintf('采用并行程序，规则数为：%d ',N)
    lengXtrain=length(Xtrain(:,1));
    R=zeros(1,lengXtrain);
    % 第1次计算2型模糊系统的输出
    tic
    parfor i11=1:lengXtrain
           [R1,R2,R(i11)]=tsk_type2(Xtrain(i11,:),M1,M2,sigma,C,S); %R(i)=1/2(R1+R2)
    end
    time11=toc
    parallel_interior_time=parallel_interior_time+time11;

    RMSE = sqrt(sum((R'-Ytrain).^2)/length(R))
    % 使用并行算法训练模型
    % 设置参数的储存变量
    M1_save=zeros(N,n,Tmax);
    M2_save=zeros(N,n,Tmax);
    C_save=zeros(N,n+1,Tmax);
    S_save=zeros(N,n+1,Tmax);
    sigma_save=zeros(N,n,Tmax);
    RMSE_parallel_save=zeros(1,Tmax);
    
    row_num=length(Xtrain(:,1));%获取当前训练数据对的组数
    if mod(row_num,core_num)==0
        lack_num=0;
    else      
        lack_num=core_num-mod(row_num,core_num);%计算填充数
        for k1=1:lack_num
            Xtrain=[Xtrain' Xtrain(k1,:)']';
            Ytrain=[Ytrain' Ytrain(k1,:)']';
        end

    end
    %R_save=zeros(Tmax,length(Ytrain)+lack_num);
    total_num=row_num+lack_num;
    circulation_num=total_num/core_num; %分配每个核心上的训练数据规模
    XX_train=zeros(core_num,length(Xtrain));
    
    time21=0;
    time22=0;
    T2=0;
    %% 使用Adam优化方法之前的参数设置，内存分配
    beta1 = 0.9;beta2=0.999;epsilon=10e-6;
    m_g_sigma = zeros(size(sigma));
    m_g_M1 = zeros(size(M1));
    m_g_M2 = zeros(size(M2));
    m_g_C = zeros(size(C));
    m_g_S = zeros(size(S));
    
    v_g_sigma = zeros(size(sigma));
    v_g_M1 = zeros(size(M1));
    v_g_M2 = zeros(size(M2));
    v_g_C = zeros(size(C));
    v_g_S = zeros(size(S)); 

    %%     
    while  T2 <Tmax %&& RMSE >=setRMSE_accept
        
           sigma_help=zeros(N,n,core_num);
           M1_help=zeros(N,n,core_num);
           M2_help=zeros(N,n,core_num);
           C_help=zeros(N,1+n,core_num);
           S_help=zeros(N,1+n,core_num);  
%            sigma22=sigma;M22=M2;M12=M1;
%            c12=C;c22=S;   
           XX_train=zeros(circulation_num,n,core_num);
           YY_train=zeros(circulation_num,1,core_num); 
           
          %% 训练数据随机重排
           X_copy=zeros(size(Xtrain));
           Y_copy=zeros(size(Ytrain));
           lengthX=length(Xtrain(:,1));
           RandID=randperm(lengthX);     
           for ji = 1:lengthX
               X_copy(ji,:) = Xtrain(RandID(ji),:);
               Y_copy(ji,:) = Ytrain(RandID(ji),:);
           end
           %% 
           tic
           parfor i22=1:core_num
                  XX_train(:,:,i22)=X_copy(((i22-1)*circulation_num+1):(i22*circulation_num),:);
                  YY_train(:,:,i22)=Y_copy(((i22-1)*circulation_num+1):(i22*circulation_num));
                  [M1_help(:,:,i22),M2_help(:,:,i22),C_help(:,:,i22),S_help(:,:,i22),sigma_help(:,:,i22)]...
                            =train_tsk_type2(XX_train(:,:,i22),YY_train(:,:,i22),M1,M2,sigma,C,S,1.5*alpha*core_num);%2
           end
           ttime21=toc;
           time21=time21+ttime21;
           
           for ii=2:core_num
               M1_help(:,:,1)=M1_help(:,:,1)+M1_help(:,:,ii);
               M2_help(:,:,1)=M2_help(:,:,1)+M2_help(:,:,ii);
               sigma_help(:,:,1)=sigma_help(:,:,1)+sigma_help(:,:,ii);
               C_help(:,:,1)=C_help(:,:,1)+C_help(:,:,ii);
               S_help(:,:,1)=S_help(:,:,1)+S_help(:,:,ii);
           end     
          
           g_sigma=sigma_help(:,:,1)/core_num;
           g_M1=M1_help(:,:,1)/core_num;
           g_M2=M2_help(:,:,1)/core_num;
           g_C=C_help(:,:,1)/core_num;
           g_S=S_help(:,:,1)/core_num;  
                 
          %% 实现Adam算法
           m_g_sigma = beta1 * m_g_sigma + (1-beta1)* g_sigma;
           m_g_M1 = beta1 * m_g_M1 + (1-beta1) * g_M1;
           m_g_M2 = beta1 * m_g_M2 + (1-beta1) * g_M2;
           m_g_C = beta1 * m_g_C + (1-beta1) * g_C;
           m_g_S = beta1 * m_g_S + (1-beta1) * g_S;
   
           v_g_sigma = beta2 * v_g_sigma + (1-beta2) * (g_sigma.*g_sigma);
           v_g_M1 = beta2 * v_g_M1 + (1-beta2) * (g_M1.*g_M1);
           v_g_M2 = beta2 * v_g_M2 + (1-beta2) * (g_M2.*g_M2);
           v_g_C = beta2 * v_g_C + (1-beta2) * (g_C.*g_C);
           v_g_S = beta2 * v_g_S + (1-beta2) * (g_S.*g_S);              
           
           m_g_sigma = m_g_sigma/(1-beta1^(T2+1));
           m_g_M1 = m_g_M1/(1-beta1^(T2+1));
           m_g_M2 = m_g_M2/(1-beta1^(T2+1));
           m_g_C = m_g_C/(1-beta1^(T2+1));
           m_g_S = m_g_S/(1-beta1^(T2+1));
           
           v_g_sigma = v_g_sigma/(1-beta2^(T2+1));
           v_g_M1 = v_g_M1/(1-beta2^(T2+1));
           v_g_M2 = v_g_M2/(1-beta2^(T2+1));
           v_g_C = v_g_C/(1-beta2^(T2+1));
           v_g_S = v_g_S/(1-beta2^(T2+1));          
           alpha1=0.02;%0.02
           sigma = g_sigma - alpha1 * (m_g_sigma./(sqrt(v_g_sigma)+epsilon));
           M1 = g_M1 - alpha1 * (m_g_M1./(sqrt(v_g_M1)+epsilon));
           M2 = g_M2 - alpha1 * (m_g_M2./(sqrt(v_g_M2)+epsilon));
           C = g_C - alpha1 * (m_g_C./(sqrt(v_g_C)+epsilon));
           S = g_S - alpha1 * (m_g_S./(sqrt(v_g_S)+epsilon));
           
          %%           
           M1_save(:,:,T2+1)=M1;
           M2_save(:,:,T2+1)=M2;
           C_save(:,:,T2+1)=C;
           S_save(:,:,T2+1)=S;
           sigma_save(:,:,T2+1)=sigma;
      
           tic
           parfor i33=1:length(Xtrain(:,1))
                  [R1,R2,R(i33)]=tsk_type2(Xtrain(i33,:),M1,M2,sigma,C,S);
           end
           ttime22=toc;
           time22=time22+ttime22;
     
           RMSE = sqrt(sum((R'-Ytrain).^2)/length(R))
           if isnan(RMSE) == 1
               RMSE = 5*setRMSE_accept
               M1 = rand(N,n);
               M2 = M1 + 1;
               sigma = rand(N,n);
               C = rand(N,1+n);
               S = rand(N,1+n); 
           end
           RMSE_parallel_save(T2+1)=RMSE;
      
           T2 = T2+1
           
%           figure(5)
%           plot(1:length(R),[R;Ytrain'],'LineWidth',1.0)
%           title(['核心数为',num2str(core_num),'规则规模为',num2str(rule_num),'，并行方法第',...
%                    num2str(T2),'次训练实际值与预测值图像'])
%           %xticks([data_size*(i-1)+1:fix(length(R)/5):data_size*(i-1)+length(R)]) 
%           xlabel('t')
%           ylabel('g(t)')
%           legend('预测值','实际值')

    end
    parallel_interior_time=parallel_interior_time+time21+time22;   
    
    wz=find(RMSE_parallel_save==min(RMSE_parallel_save(RMSE_parallel_save~=0)));
    wz=wz(1);
    Rs=zeros(1,length(Xtrain(:,1)));
    
    tic
    parfor i=1:length(Xtrain(:,1))
            [R11,R22,Rs(i)]=tsk_type2(Xtrain(i,:),M1_save(:,:,wz),M2_save(:,:,wz),sigma_save(:,:,wz),C_save(:,:,wz),S_save(:,:,wz));
    end
    time12=toc;
    parallel_interior_time=parallel_interior_time+time12;
    
    RMSE = sqrt(sum((Rs'-Ytrain).^2)/length(Rs))
    
    end_parallel_time=clock;
    parallel_total_time=etime(end_parallel_time,star_parallel_time)  
    %%
    SAVE_parallel_M1(core_num,1:length(M1_save(:,1,wz)),:,rule_num)    = M1_save(:,:,wz);
    SAVE_parallel_M2(core_num,1:length(M2_save(:,1,wz)),:,rule_num)    = M2_save(:,:,wz);
    SAVE_parallel_sigma(core_num,1:length(sigma_save(:,1,wz)),:,rule_num) = sigma_save(:,:,wz);
    SAVE_parallel_C(core_num,1:length(C_save(:,1,wz)),:,rule_num)      = C_save(:,:,wz);
    SAVE_parallel_S(core_num,1:length(S_save(:,1,wz)),:,rule_num)      = S_save(:,:,wz);
    SAVE_parallel_Ytrain(core_num,1:length(Ytrain),DATA_SIZE,rule_num)       = Ytrain;
    SAVE_parallel_Ypridect(core_num,1:length(Rs),DATA_SIZE,rule_num) = Rs;   

    %% 记录并行总时间，并行内部时间，串行时间，计算加速比、并行效率,串行比例
    SAVEparallel_total_time(core_num,i,rule_num)    = parallel_total_time;      
    SAVEparallel_interior_time(core_num,i,rule_num) = parallel_interior_time;         
    SAVEspeedUP_ratio(core_num,i,rule_num)          = SAVEserial_total_time(rule_num,i)/parallel_total_time;
    SAVEparallel_efficiency(core_num,i,rule_num)    = SAVEserial_total_time(rule_num,i)/(core_num*parallel_total_time);
    SAVEserial_proportion(core_num,i,rule_num)     = (parallel_total_time - parallel_interior_time)/parallel_total_time;   
    SAVE_RMSE_parallel(core_num,:,i,rule_num) = RMSE_parallel_save;
    
end       
    delete(gcp('nocreate')) %每一轮循环后关闭MATLAB并行计算池，下一轮开启不同核心数的并行计算池 
end  

%% 保存所有需要的数据到SAVE_DATA.mat文件中，以便于列表或做图
save('A','SAVEparallel_total_time','SAVEparallel_interior_time',...
    'SAVEserial_total_time','SAVEspeedUP_ratio','SAVEparallel_efficiency',...
    'SAVE_RMSE_serial','SAVE_RMSE_parallel','SAVEserial_proportion')
save('B','SAVE_parallel_M1','SAVE_parallel_M2','SAVE_parallel_sigma',...
    'SAVE_parallel_C','SAVE_parallel_S','SAVE_parallel_Ytrain','SAVE_parallel_Ypridect')
save('C','SAVE_serial_M1','SAVE_serial_M2','SAVE_serial_sigma',...
    'SAVE_serial_C','SAVE_serial_S','SAVE_serial_Ytest','SAVE_serial_Ypridect')



