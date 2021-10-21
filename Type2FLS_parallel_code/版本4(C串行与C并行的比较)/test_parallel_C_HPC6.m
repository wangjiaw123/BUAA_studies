close all
clear
clc

%% ��������
%star_parallel_time=clock;
tao=31;
Dt=tao;
N=40020;%���õ���N��Mackey_GlassĬ��ʱ϶���Ϊ0.1;
[y,t]=Mackey_Glass(N,tao);%����-����(Runge-Kutta)�������Mackey Glass����
%% ԭʼͼ��
% figure(1)
% subplot(1,2,1)% ����ֱ�� 
% plot(t,y,'LineWidth',1.0);
% xlabel('t')
% ylabel('s(t)')
% title('����ֱ�� ');
% 
% subplot(1,2,2)% ��ͼʱ�� 
% plot(y(Dt*10+10000+1:end),y(10001:end-Dt*10),'LineWidth',0.5)
% xlabel('s(t)')
% ylabel('s(t-\tau)')
% title('��ͼʱ��');

n_train = 4000;
x_star = zeros(1,n_train);
for i = 1:n_train
    x_star(i) = y(1+i*10);
end

%N = 40;%���ù�����Ŀ20
%% 

%% ȫ�ֱ���
n = 4;%ģ������ǰ������Ϊ4
alpha =0.4;%KMʱȡ0.4
setRMSE_accept=0.04;  %���ÿɽ��ܵ����
%% �˴�������Ҫ�������ֵ���󣬶�ά����һ��ʼ��С��Ҫ�������!
max_select_train = 2100;
data_size = 500;%1000
RULE=[20,40];%[20,60]
RULE_NUMBER = length(RULE);%10
RULE_MAX = max(RULE);
Monte_carlo=1;%���ؿ���ʵ��20
CORE_NUMBER =6;%2:2:24      
%CORE_NUMBER_WIDTH=1;
DATA_SIZE = 1;%7  %max_select_train/data_size=7                              
Tmax = 10;  %40                                       %����ģ��ѵ���е�����������

                                
SAVEparallel_total_time_C    = zeros(CORE_NUMBER,DATA_SIZE,RULE_NUMBER,Monte_carlo);%1ά�Ƿ���ĺ������ĸ�����2ά�����ݹ�ģ���ĸ�����3ά�ǹ������ĸ���             
SAVEparallel_interior_time_C = zeros(CORE_NUMBER,DATA_SIZE,RULE_NUMBER,Monte_carlo);%�ֱ��ǲ�����ʱ�䣬�����ڲ�ʱ�䣬����ʱ�䣬������ٱȡ�����Ч�ʡ����б�������ĳ�ʼ��        
SAVE_RMSE_parallel_C = zeros(CORE_NUMBER,Tmax,DATA_SIZE,RULE_NUMBER,Monte_carlo);   %1ά�Ƿ���ĺ������ĸ�����2ά��Tmax(��һά���ڴ洢RMSE),3ά�����ݹ�ģ���ĸ�����4ά�ǹ������ĸ���

SAVE_parallel_C_M1    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_MAX,n,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_parallel_C_M2    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_MAX,n,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_parallel_C_sigma = zeros(CORE_NUMBER,RULE_NUMBER*RULE_MAX,n,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_parallel_C_c1    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_MAX,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_parallel_C_c2    = zeros(CORE_NUMBER,RULE_NUMBER*RULE_MAX,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_parallel_C_Ytest = zeros(CORE_NUMBER,max_select_train,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_parallel_Ypridect_C = zeros(CORE_NUMBER,data_size,DATA_SIZE,RULE_NUMBER,Monte_carlo);

SAVEparallel_total_time_C_MEAN    = zeros(CORE_NUMBER,DATA_SIZE,RULE_NUMBER);      
SAVEparallel_interior_time_C_MEAN = zeros(CORE_NUMBER,DATA_SIZE,RULE_NUMBER);        
SAVE_RMSE_parallel_C_MEAN           = zeros(CORE_NUMBER,Tmax,DATA_SIZE,RULE_NUMBER);
SAVE_parallel_Ypridect_C_MEAN       = zeros(CORE_NUMBER,data_size,DATA_SIZE,RULE_NUMBER);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SAVEserial_total_time_C = zeros(RULE_NUMBER,DATA_SIZE,Monte_carlo);        %1ά�ǹ������ĸ�����2ά�����ݹ�ģ���ĸ���
SAVE_RMSE_serial_C    = zeros(RULE_NUMBER,Tmax,DATA_SIZE,Monte_carlo);     %ά�ǹ������ĸ�����2ά��Tmax(��һά���ڴ洢RMSE),3ά�����ݹ�ģ���ĸ���

SAVE_serial_M1_C    = zeros(RULE_NUMBER*RULE_MAX,n,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_serial_M2_C    = zeros(RULE_NUMBER*RULE_MAX,n,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_serial_sigma_C = zeros(RULE_NUMBER*RULE_MAX,n,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_serial_c1_C    = zeros(RULE_NUMBER*RULE_MAX,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_serial_c2_C    = zeros(RULE_NUMBER*RULE_MAX,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_serial_Ytest_C = zeros(max_select_train,DATA_SIZE,RULE_NUMBER,Monte_carlo);
SAVE_serial_Ypridect_C = zeros(data_size,DATA_SIZE,RULE_NUMBER,Monte_carlo);         


SAVEserial_total_time_C_MEAN = zeros(RULE_NUMBER,DATA_SIZE);
SAVE_RMSE_serial_C_MEAN    = zeros(RULE_NUMBER,Tmax,DATA_SIZE);  
SAVE_serial_Ypridect_C_MEAN = zeros(data_size,DATA_SIZE,RULE_NUMBER);   


%% ���г��򲿷�

for MTC=1:Monte_carlo

for rule_num=1:RULE_NUMBER
    
    N=RULE(rule_num);          %������Ŀ��ѭ��
    
for i=1:DATA_SIZE                   %�����ģ��ѭ��
   %% ʱ���¼������ʼ��
    serial_total_time = 0;   
    %% ѡ��ѵ���������������
    % ֻ��һ����ģ������ѵ��ģ���߼�ϵͳ����Ԥ���500�����������磺��1��500ѵ����
    % Ԥ��501��1000��...,��1��3500ѵ����Ԥ��3501��4000
    X_train = [x_star(1:i*data_size-4);x_star(2:i*data_size-3);x_star(3:i*data_size-2);x_star(4:i*data_size-1);]';
    Y_train = x_star(5:i*data_size)';
    X_test = [x_star(i*data_size-3:(i+1)*data_size-4);x_star(i*data_size-2:(i+1)*data_size-3);...
            x_star(i*data_size-1:(i+1)*data_size-2);x_star(i*data_size:(i+1)*data_size-1)]';
    Y_test = x_star(i*data_size+1:(i+1)*data_size)';
    
    %% ѵ��������ʼ��
    M1 = rand(N,n);
    M2 = M1 + 0.8;
    sigma = rand(N,n);
    c1 = rand(N,1);
    c2 = c1+1;    
    
    %% ���з���
    fprintf('������Ϊ��%d,��ʼ�� %d �ֹ�ģ��ѵ��(���ô��г���)��ѵ�����ݹ�ģΪ:%d ,Ԥ�����ݹ�ģΪ:%d \n',N,i,data_size*i,data_size);
    star_serial_time=clock;
    R=[];
    for i1=1:length(X_train(:,1))
          %[R1,R2,R(i1)]=sfls_type2(X_train(i1,:),M1,M2,sigma,c1,c2); %R(i)=1/2(R1+R2)
          [R1,R2,R(i1)]=sfls_type2_C(X_train(i1,:),M1,M2,c1,c2,sigma); %R(i)=1/2(R1+R2)
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
           %[M1,M2,c1,c2,sigma,I2l,I2u,I1u,I1l]=train_sfls_type2(X_train,Y_train,M1,M2,sigma,c1,c2,alpha);
           [M1,M2,c1,c2,sigma]=Train_sfls_type2_C(X_train,Y_train,M1,M2,c1,c2,sigma,(1-T1/Tmax+0.005)*alpha);%
           for i2=1:length(X_train(:,1))
                  %[R1,R2,R(i2)]=sfls_type2(X_train(i2,:),M1,M2,sigma,c1,c2);
                  [R1,R2,R(i2)]=sfls_type2_C(X_train(i2,:),M1,M2,c1,c2,sigma);
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
%            title(['�����ģΪ',num2str(rule_num),',�����ģ��Ϊ',num2str(i),...
%                '�����з�����',num2str(T1),'��ѵ��ʵ��ֵ��Ԥ��ֵͼ��','(\tau =',num2str(tao),')'])
%            %xticks([data_size*(i-1)+1:fix(length(R)/5):data_size*(i-1)+length(R)]) 
%            xlabel('t')
%            ylabel('s(t)')
%            legend('Ԥ��ֵ','ʵ��ֵ')      
                        
    end    

    wz_serial=find(RMSE_serial_save==min(RMSE_serial_save(RMSE_serial_save~=0)));
    wz_serial=wz_serial(1);
    Rs_serial=zeros(1,length(X_test(:,1)));   
    for ik=1:length(X_test(:,1))
%              [R11,R22,Rs_serial(ik)]=sfls_type2(X_test(ik,:),M1_serial_save(:,:,wz_serial),...
%                  M2_serial_save(:,:,wz_serial),sigma_serial_save(:,:,wz_serial),...
%                  c1_serial_save(:,:,wz_serial),c2_serial_save(:,:,wz_serial));
           [R11,R22,Rs_serial(ik)]=sfls_type2_C(X_test(ik,:),M1_serial_save(:,:,wz_serial),...
               M2_serial_save(:,:,wz_serial),c1_serial_save(:,:,wz_serial)...
               ,c2_serial_save(:,:,wz_serial),sigma_serial_save(:,:,wz_serial));
    end
    RMSE = sqrt(sum((Rs_serial'-Y_test).^2)/length(Rs_serial))
    
    end_serial_time=clock;
    serial_total_time=etime(end_serial_time,star_serial_time)  
    SAVEserial_total_time_C(rule_num,i,MTC)      = serial_total_time;   
    
    SAVE_serial_M1_C(1:length(M1_serial_save(:,1,wz_serial)),:,i,rule_num,MTC)    = M1_serial_save(:,:,wz_serial);
    SAVE_serial_M2_C(1:length(M2_serial_save(:,1,wz_serial)),:,i,rule_num,MTC)    = M2_serial_save(:,:,wz_serial);
    SAVE_serial_sigma_C(1:length(sigma_serial_save(:,1,wz_serial)),:,i,rule_num,MTC) = sigma_serial_save(:,:,wz_serial);
    SAVE_serial_c1_C(1:length(c1_serial_save(:,1,wz_serial)),i,rule_num,MTC)      = c1_serial_save(:,:,wz_serial);
    SAVE_serial_c2_C(1:length(c2_serial_save(:,1,wz_serial)),i,rule_num,MTC)      = c2_serial_save(:,:,wz_serial);
    SAVE_serial_Ytest_C(1:length(Y_test),i,rule_num,MTC)       = Y_test;
    SAVE_serial_Ypridect_C(1:length(Rs_serial),i,rule_num,MTC) = Rs_serial;
    SAVE_RMSE_serial_C(rule_num,:,i,MTC)    = RMSE_serial_save; 
    
%     figure(3)
%     subplot(1,2,1)
%     plot(i*data_size+1:i*data_size+length(Rs_serial),Rs_serial,'LineWidth',1.0)
%     xlabel('t')
%     ylabel('s(t)')
%     title('����ֱ�� ');
%     subplot(1,2,2)
%     plot(Rs_serial(tao+1:end),Rs_serial(1:end-tao),'LineWidth',0.5)
%     xlabel('s(t)')
%     ylabel('s(t-\tau)')
%     title('��ͼʱ��');
% 
%     figure(4)
%     plot(i*data_size+1+length(R):i*data_size+length(Rs_serial)+length(R),[Rs_serial;Y_test'],'LineWidth',1.0)
%     xlabel('t')
%     ylabel('s(t)')
%     title(['���з���ʵ��ֵ������Ԥ��ֵͼ��','(\tau =',num2str(tao),')'])
%     legend('Ԥ��ֵ','ʵ��ֵ')    
  
end
        
end


end    %end Monte_carlo

for MTC1=1:Monte_carlo
    SAVEserial_total_time_C_MEAN =SAVEserial_total_time_C_MEAN + SAVEserial_total_time_C(:,:,MTC1);
    SAVE_RMSE_serial_C_MEAN    = SAVE_RMSE_serial_C_MEAN + SAVE_RMSE_serial_C(:,:,:,MTC1);  
    SAVE_serial_Ypridect_C_MEAN =SAVE_serial_Ypridect_C_MEAN + SAVE_serial_Ypridect_C(:,:,:,MTC1);   
end

SAVEserial_total_time_C_MEAN =SAVEserial_total_time_C_MEAN/Monte_carlo;
SAVE_RMSE_serial_C_MEAN    = SAVE_RMSE_serial_C_MEAN/Monte_carlo;
SAVE_serial_Ypridect_C_MEAN =SAVE_serial_Ypridect_C_MEAN/Monte_carlo;






%% ���г��򲿷�
for MC=1:Monte_carlo


for core_num=2:2:CORE_NUMBER %   %��������ѭ��
    
    parpool(core_num);          %��MATLAB�п���������Ϊcore_num�Ĳ��м����  
    poolobj = gcp('nocreate')
%    addAttachedFiles(poolobj,{'/gs/home/sy1909131/Train_sfls_type2_C.mexa64','/gs/home/sy1909131/sfls_type2_C.mexa64'})
%     addAttachedFiles(poolobj,{'C:\Users\wangjiawen\Desktop\�½��ļ��� (2)\Train_sfls_type2_C.mexw64',...
%                  'C:\Users\wangjiawen\Desktop\�½��ļ��� (2)\sfls_type2_C.mexw64'})
    addAttachedFiles(poolobj,{'E:\������\��Ҫ�ļ�\������ʦ���ļ�\����ģ���߼�ϵͳ�Ĳ���ʵ��\2��ģ���߼�ϵͳ---���г���\����(��ȷ��)\�汾4(C������C���еıȽ�)\Train_sfls_type2_C.mexw64',...
                'E:\������\��Ҫ�ļ�\������ʦ���ļ�\����ģ���߼�ϵͳ�Ĳ���ʵ��\2��ģ���߼�ϵͳ---���г���\����(��ȷ��)\�汾4(C������C���еıȽ�)\sfls_type2_C.mexw64'})   
            
for rule_num=1:RULE_NUMBER
    
    N=RULE(rule_num);                     %������Ŀ��ѭ��
    %N=RULE(rule_num); 
for i=1:DATA_SIZE                      %�����ģ��ѭ��
   %% ʱ���¼������ʼ��
    parallel_total_time = 0;
    parallel_interior_time = 0;
    %% ѡ��ѵ���������������
    % ֻ��һ����ģ������ѵ��ģ���߼�ϵͳ����Ԥ���500�����������磺��1��500ѵ����
    % Ԥ��501��1000��...,��3001��3500ѵ����Ԥ��3501��4000    X_train = [x_star(1:i*data_size-4);x_star(2:i*data_size-3);x_star(3:i*data_size-2);x_star(4:i*data_size-1);]';
    X_train = [x_star(1:i*data_size-4);x_star(2:i*data_size-3);x_star(3:i*data_size-2);x_star(4:i*data_size-1);]';
    Y_train = x_star(5:i*data_size)';
    X_test = [x_star(i*data_size-3:(i+1)*data_size-4);x_star(i*data_size-2:(i+1)*data_size-3);...
            x_star(i*data_size-1:(i+1)*data_size-2);x_star(i*data_size:(i+1)*data_size-1)]';
    Y_test = x_star(i*data_size+1:(i+1)*data_size)';
     
    %% ѵ��������ʼ��
    M1 = rand(N,n);
    M2 = M1 + 1;%0.3
    sigma = rand(N,n);
    c1 = rand(N,1);
    c2 = c1+1.5; %0.55
    %% ���з���
    star_parallel_time=clock;%ֱ��������������㲢�г����ܵ�����ʱ��
    fprintf('������Ϊ��%d ,��ʼ�� %d �ֹ�ģ��ѵ��(���ò��г���)��ѵ�����ݹ�ģΪ:%d ,Ԥ�����ݹ�ģΪ:%d \n',N,i,data_size*i,data_size)
    R=[];
    % ��1�μ���2��ģ��ϵͳ�����
    tic
    parfor i11=1:length(X_train(:,1))
    %for i11=1:length(X_train(:,1))
           %[R1,R2,R(i11)]=sfls_type2(X_train(i11,:),M1,M2,sigma,c1,c2); %R(i)=1/2(R1+R2)
           [R1,R2,R(i11)]=sfls_type2_C(X_train(i11,:),M1,M2,c1,c2,sigma);
    end
    time11=toc;
    parallel_interior_time=parallel_interior_time+time11;

    RMSE = sqrt(sum((R'-Y_train).^2)/length(R))
    % ʹ�ò����㷨ѵ��ģ��
    % ���ò����Ĵ������
    M1_save=zeros(N,n,Tmax);
    M2_save=zeros(N,n,Tmax);
    c1_save=zeros(N,1,Tmax);
    c2_save=zeros(N,1,Tmax);
    sigma_save=zeros(N,n,Tmax);
    RMSE_parallel_save=zeros(1,Tmax);
    
    row_num=length(X_train(:,1));%��ȡ��ǰѵ�����ݶԵ�����
    lack_num=core_num-mod(row_num,core_num);%���������
    R_save=zeros(Tmax,length(Y_train)+lack_num);
    for k=1:lack_num
        X_train=[X_train' X_train(k,:)']';
        Y_train=[Y_train' Y_train(k,:)']';
    end
    total_num=row_num+lack_num;
    circulation_num=total_num/core_num; %����ÿ�������ϵ�ѵ�����ݹ�ģ
    XX_train=zeros(core_num,length(X_train));
    
    time21=0;
    time22=0;
    T2=0;
    %% ʹ��Adam�Ż�����֮ǰ�Ĳ������ã��ڴ����
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
           
          %% ѵ�������������
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
           for i22=1:core_num
                  XX_train(:,:,i22)=X_copy(((i22-1)*circulation_num+1):(i22*circulation_num),:);
                  YY_train(:,:,i22)=Y_copy(((i22-1)*circulation_num+1):(i22*circulation_num));
%                   [M1_help(:,:,i22),M2_help(:,:,i22),c1_help(:,:,i22),c2_help(:,:,i22),sigma_help(:,:,i22)]...
%                             =train_sfls_type2(XX_train(:,:,i22),YY_train(:,:,i22),M1,M2,sigma,c1,c2,alpha*core_num);%2
                  [M1_help(:,:,i22),M2_help(:,:,i22),c1_help(:,:,i22),c2_help(:,:,i22),sigma_help(:,:,i22)]...      
                       =Train_sfls_type2_C(XX_train(:,:,i22),YY_train(:,:,i22),M1,M2,c1,c2,sigma,2*alpha);%����C����ѧϰ�ʲ��ܵ�̫��1*alpha����
           
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
          %% ʵ��Adam�㷨
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
           %for i33=1:length(X_train(:,1))
                  %[R1,R2,R(i33)]=sfls_type2(X_train(i33,:),M1,M2,sigma,c1,c2);
                  [R1,R2,R(i33)]=sfls_type2_C(X_train(i33,:),M1,M2,c1,c2,sigma);
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
%            title(['������Ϊ',num2str(core_num),'�����ģΪ',num2str(rule_num),',�����ģ��Ϊ',num2str(i),...
%               '�����з�����',num2str(T2),'��ѵ��ʵ��ֵ��Ԥ��ֵͼ��','(\tau =',num2str(tao),')'])
%            %xticks([data_size*(i-1)+1:fix(length(R)/5):data_size*(i-1)+length(R)]) 
%            xlabel('t')
%            ylabel('s(t)')
%            legend('Ԥ��ֵ','ʵ��ֵ')

    end
    parallel_interior_time=parallel_interior_time+time21+time22;   
    
    wz=find(RMSE_parallel_save==min(RMSE_parallel_save(RMSE_parallel_save~=0)));
    wz=wz(1);
    Rs=zeros(1,length(X_test(:,1)));
    
    tic
    parfor i=1:length(X_test(:,1))
    %for i=1:length(X_test(:,1))
            %[R11,R22,Rs(i)]=sfls_type2(X_test(i,:),M1_save(:,:,wz),M2_save(:,:,wz),sigma_save(:,:,wz),c1_save(:,:,wz),c2_save(:,:,wz));
            [R11,R22,Rs(i)]=sfls_type2_C(X_test(i,:),M1_save(:,:,wz),M2_save(:,:,wz),c1_save(:,:,wz),c2_save(:,:,wz),sigma_save(:,:,wz));
    end
    time12=toc;
    parallel_interior_time=parallel_interior_time+time12;
    
    RMSE = sqrt(sum((Rs'-Y_test).^2)/length(Rs))
    
    end_parallel_time=clock;
    parallel_total_time=etime(end_parallel_time,star_parallel_time)  
    
    SAVE_parallel_C_M1(core_num,1:length(M1_save(:,1,wz)),:,i,rule_num,MC)    = M1_save(:,:,wz);
    SAVE_parallel_C_M2(core_num,1:length(M2_save(:,1,wz)),:,i,rule_num,MC)    = M2_save(:,:,wz);
    SAVE_parallel_C_sigma(core_num,1:length(sigma_save(:,1,wz)),:,i,rule_num,MC) = sigma_save(:,:,wz);
    SAVE_parallel_C_c1(core_num,1:length(c1_save(:,1,wz)),i,rule_num,MC)      = c1_save(:,:,wz);
    SAVE_parallel_C_c2(core_num,1:length(c2_save(:,1,wz)),i,rule_num,MC)      = c2_save(:,:,wz);
    SAVE_parallel_C_Ytest(core_num,1:length(Y_test),i,rule_num,MC)       = Y_test;
    SAVE_parallel_Ypridect_C(core_num,1:length(Rs),i,rule_num,MC) = Rs;   

%     figure(6)
%     subplot(1,2,1)
%     plot(i*data_size+1:i*data_size+length(Rs),Rs,'LineWidth',1.0)
%     xlabel('t')
%     ylabel('s(t)')
%     title('����ֱ�� ');
%     subplot(1,2,2)
%     plot(Rs(tao+1:end),Rs(1:end-tao),'LineWidth',0.5)
%     xlabel('s(t)')
%     ylabel('s(t-\tau)')
%     title('��ͼʱ��');
% 
%     figure(7)
%     plot(i*data_size+1+length(R):i*data_size+length(Rs)+length(R),[Rs;Y_test'],'LineWidth',1.0)
%     xlabel('t')
%     ylabel('s(t)')
%     title(['���з���ʵ��ֵ������Ԥ��ֵͼ��','(\tau =',num2str(tao),')'])
%     legend('Ԥ��ֵ','ʵ��ֵ')

    %% ��¼������ʱ�䣬�����ڲ�ʱ�䣬����ʱ�䣬������ٱȡ�����Ч��,���б���
    SAVEparallel_total_time_C(core_num,i,rule_num,MC)    = parallel_total_time;      
    SAVEparallel_interior_time_C(core_num,i,rule_num,MC) = parallel_interior_time;          
    SAVE_RMSE_parallel_C(core_num,:,i,rule_num,MC) = RMSE_parallel_save;
end       
    
end  
   delete(gcp('nocreate')) %ÿһ��ѭ����ر�MATLAB���м���أ���һ�ֿ�����ͬ�������Ĳ��м����
end

end

for MCMC=1:Monte_carlo
    SAVEparallel_total_time_C_MEAN    = SAVEparallel_total_time_C_MEAN+SAVEparallel_total_time_C(:,:,:,MCMC);   
    SAVEparallel_interior_time_C_MEAN = SAVEparallel_interior_time_C_MEAN+SAVEparallel_interior_time_C(:,:,:,MCMC);        
    SAVE_RMSE_parallel_C_MEAN         = SAVE_RMSE_parallel_C_MEAN+SAVE_RMSE_parallel_C(:,:,:,:,MCMC);
    SAVE_parallel_Ypridect_C_MEAN       =SAVE_parallel_Ypridect_C_MEAN +SAVE_parallel_Ypridect_C(:,:,:,:,MCMC);
end
SAVEparallel_total_time_C_MEAN    =SAVEparallel_total_time_C_MEAN/Monte_carlo;   
SAVEparallel_interior_time_C_MEAN =SAVEparallel_interior_time_C_MEAN/Monte_carlo ;        
SAVE_RMSE_parallel_C_MEAN         =SAVE_RMSE_parallel_C_MEAN/Monte_carlo;
SAVE_parallel_Ypridect_C_MEAN       =SAVE_parallel_Ypridect_C_MEAN/Monte_carlo;


%% ����������Ҫ�����ݵ�SAVE_DATA.mat�ļ��У��Ա����б����ͼ

save('C_SAVEdata_parallel_parameter','SAVEparallel_total_time_C','SAVEparallel_interior_time_C',...
    'SAVE_RMSE_parallel_C','SAVE_parallel_C_M1','SAVE_parallel_C_M2','SAVE_parallel_C_sigma',...
    'SAVE_parallel_C_c1','SAVE_parallel_C_c2','SAVE_parallel_C_Ytest','SAVE_parallel_Ypridect_C')
save('C_SAVEdata_serial_parameter','SAVEserial_total_time_C','SAVE_RMSE_serial_C','SAVE_serial_M1_C',...
    'SAVE_serial_M2_C','SAVE_serial_sigma_C','SAVE_serial_c1_C','SAVE_serial_c2_C','SAVE_serial_Ytest_C','SAVE_serial_Ypridect_C')
save('C_SAVEdata_MEAN','SAVEparallel_total_time_C_MEAN','SAVEparallel_interior_time_C_MEAN',...
    'SAVE_RMSE_parallel_C_MEAN','SAVE_parallel_Ypridect_C_MEAN','SAVEserial_total_time_C_MEAN',...
    'SAVE_RMSE_serial_C_MEAN','SAVE_serial_Ypridect_C_MEAN')



