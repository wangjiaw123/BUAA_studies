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
RULE=[20,40,60,80];%[20,60]
RULE_NUMBER = length(RULE);%10
RULE_MAX = max(RULE);
Monte_carlo=10;%���ؿ���ʵ��20
CORE_NUMBER =24;%2:2:24      
%CORE_NUMBER_WIDTH=1;
DATA_SIZE = 7;%7  %max_select_train/data_size=7                              
Tmax = 30;  %40                                       %����ģ��ѵ���е�����������

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
    fprintf('���ؿ��޴���Ϊ��%d,������Ϊ��%d,��ʼ�� %d �ֹ�ģ��ѵ��(���ô��г���)��ѵ�����ݹ�ģΪ:%d ,Ԥ�����ݹ�ģΪ:%d \n',MTC,N,i,data_size*i,data_size);
    star_serial_time=clock;
    R=[];
    for i1=1:length(X_train(:,1))
          [R1,R2,R(i1)]=sfls_type2(X_train(i1,:),M1,M2,sigma,c1,c2); %R(i)=1/2(R1+R2)
          %[R1,R2,R(i1)]=sfls_type2_C(X_train(i1,:),M1,M2,c1,c2,sigma); %R(i)=1/2(R1+R2)
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
           [M1,M2,c1,c2,sigma,I2l,I2u,I1u,I1l]=train_sfls_type2(X_train,Y_train,M1,M2,sigma,c1,c2,alpha);
           %[M1,M2,c1,c2,sigma]=Train_sfls_type2_C(X_train,Y_train,M1,M2,c1,c2,sigma,(1-T1/Tmax+0.005)*alpha);%
           for i2=1:length(X_train(:,1))
                  [R1,R2,R(i2)]=sfls_type2(X_train(i2,:),M1,M2,sigma,c1,c2);
                  %[R1,R2,R(i2)]=sfls_type2_C(X_train(i2,:),M1,M2,c1,c2,sigma);
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
             [R11,R22,Rs_serial(ik)]=sfls_type2(X_test(ik,:),M1_serial_save(:,:,wz_serial),...
                 M2_serial_save(:,:,wz_serial),sigma_serial_save(:,:,wz_serial),...
                 c1_serial_save(:,:,wz_serial),c2_serial_save(:,:,wz_serial));
%            [R11,R22,Rs_serial(ik)]=sfls_type2_C(X_test(ik,:),M1_serial_save(:,:,wz_serial),...
%                M2_serial_save(:,:,wz_serial),c1_serial_save(:,:,wz_serial)...
%                ,c2_serial_save(:,:,wz_serial),sigma_serial_save(:,:,wz_serial));
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


save('data_serial_MonteCarlo','SAVEserial_total_time_C','SAVE_RMSE_serial_C','SAVE_serial_M1_C',...
    'SAVE_serial_M2_C','SAVE_serial_sigma_C','SAVE_serial_c1_C','SAVE_serial_c2_C','SAVE_serial_Ytest_C',...
    'SAVE_serial_Ypridect_C');
save('data_serial_MEAN','SAVEserial_total_time_C_MEAN','SAVE_RMSE_serial_C_MEAN','SAVE_serial_Ypridect_C_MEAN');



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







