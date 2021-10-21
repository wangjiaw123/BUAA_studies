
close all
clear
clc

%% ��������
%star_parallel_time=clock;
tao=13;
Dt=tao;
N=40020;%���õ���N��Mackey_GlassĬ��ʱ϶���Ϊ0.1;
[y,t]=Mackey_Glass(N,tao);%����-����(Runge-Kutta)�������Mackey Glass����
%% ԭʼͼ��
figure(1)
subplot(1,2,1)% ����ֱ�� 
plot(t,y,'LineWidth',1.0);
xlabel('t')
ylabel('s(t)')
title('����ֱ�� ');

subplot(1,2,2)% ��ͼʱ�� 
plot(y(Dt*10+10000+1:end),y(10001:end-Dt*10),'LineWidth',0.5)
xlabel('s(t)')
ylabel('s(t-\tau)')
title('��ͼʱ��');

n_train = 4000;
x_star = zeros(1,n_train);
for i = 1:n_train
    x_star(i) = y(1+i*10);
end


%N = 40;%���ù�����Ŀ20
%% 
max_select_train = 3500;
data_size = 500;

%% ȫ�ֱ���
n = 4;%ģ������ǰ������Ϊ4
alpha =0.4;%KMʱȡ0.4
setRMSE_accept=0.04;  %���ÿɽ��ܵ����
%% �˴�������Ҫ�������ֵ���󣬶�ά����һ��ʼ��С��Ҫ������ˣ�
RULR_NUMBER = 10;
RULE_WIDTH = 10;
CORE_NUMBER =36;      %�鿴���в��֣���������ʹ���ǵݼ��ģ�
%CORE_NUMBER_WIDTH=1;
DATA_SIZE = 7;  %max_select_train/data_size=7                              
Tmax = 160;                                         %����ģ��ѵ���е�����������

Zero1 = zeros(CORE_NUMBER,DATA_SIZE,RULR_NUMBER);  %1ά�Ƿ���ĺ������ĸ�����2ά�����ݹ�ģ���ĸ�����3ά�ǹ������ĸ���
SAVEparallel_total_time_A    = Zero1;              %�ֱ��ǲ�����ʱ�䣬�����ڲ�ʱ�䣬����ʱ�䣬������ٱȡ�����Ч�ʡ����б�������ĳ�ʼ��
SAVEparallel_interior_time_B = Zero1;         
SAVEspeedUP_ratio_D          = Zero1;
SAVEparallel_efficiency_E    = Zero1;
SAVESserial_proportion_H     = Zero1;

SAVE_RMSE__parallel_G = zeros(CORE_NUMBER,Tmax,DATA_SIZE,RULR_NUMBER);   %1ά�Ƿ���ĺ������ĸ�����2ά��Tmax(��һά���ڴ洢RMSE),3ά�����ݹ�ģ���ĸ�����4ά�ǹ������ĸ���

SAVE_parallel_M1    = zeros(CORE_NUMBER,RULR_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULR_NUMBER);
SAVE_parallel_M2    = zeros(CORE_NUMBER,RULR_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULR_NUMBER);
SAVE_parallel_sigma = zeros(CORE_NUMBER,RULR_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULR_NUMBER);
SAVE_parallel_c1    = zeros(CORE_NUMBER,RULR_NUMBER*RULE_WIDTH,DATA_SIZE,RULR_NUMBER);
SAVE_parallel_c2    = zeros(CORE_NUMBER,RULR_NUMBER*RULE_WIDTH,DATA_SIZE,RULR_NUMBER);
SAVE_parallel_Ytest = zeros(CORE_NUMBER,max_select_train,DATA_SIZE,RULR_NUMBER);
SAVE_parallel_Ypridect = zeros(CORE_NUMBER,data_size,DATA_SIZE,RULR_NUMBER);

SAVEserial_total_time_C = zeros(RULR_NUMBER,DATA_SIZE);        %1ά�ǹ������ĸ�����2ά�����ݹ�ģ���ĸ���
SAVE_RMSE_serial_F    = zeros(RULR_NUMBER,Tmax,DATA_SIZE);     %ά�ǹ������ĸ�����2ά��Tmax(��һά���ڴ洢RMSE),3ά�����ݹ�ģ���ĸ���

SAVE_serial_M1    = zeros(RULR_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULR_NUMBER);
SAVE_serial_M2    = zeros(RULR_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULR_NUMBER);
SAVE_serial_sigma = zeros(RULR_NUMBER*RULE_WIDTH,n,DATA_SIZE,RULR_NUMBER);
SAVE_serial_c1    = zeros(RULR_NUMBER*RULE_WIDTH,DATA_SIZE,RULR_NUMBER);
SAVE_serial_c2    = zeros(RULR_NUMBER*RULE_WIDTH,DATA_SIZE,RULR_NUMBER);
SAVE_serial_Ytest = zeros(max_select_train,DATA_SIZE,RULR_NUMBER);
SAVE_serial_Ypridect = zeros(data_size,DATA_SIZE,RULR_NUMBER);         



%% ���г��򲿷�
for rule_num=1:RULR_NUMBER
    N=RULE_WIDTH*rule_num;          %������Ŀ��ѭ��

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
    M2 = M1 + 1;
    sigma = rand(N,n);
    c1 = rand(N,1);
    c2 = c1+1.5;    
    
    %% ���з���
    fprintf('������Ϊ��%d,��ʼ�� %d �ֹ�ģ��ѵ��(���ô��г���)��ѵ�����ݹ�ģΪ:%d ,Ԥ�����ݹ�ģΪ:%d \n',N,i,data_size*i,data_size);
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
    
    while  T1 <Tmax && RMSE  >=setRMSE_accept
           [M1,M2,c1,c2,sigma,I2l,I2u,I1u,I1l]=train_sfls_type2(X_train,...
               Y_train,M1,M2,sigma,c1,c2,alpha);
           
           for i2=1:length(X_train(:,1))
                  [R1,R2,R(i2)]=sfls_type2(X_train(i2,:),M1,M2,sigma,c1,c2);
           end     
           RMSE = sqrt(sum((R'-Y_train).^2)/length(R))
           RMSE_serial_save(T1+1)=RMSE; 
           
           M1_serial_save(:,:,T1+1)=M1;
           M2_serial_save(:,:,T1+1)=M2;
           sigma_serial_save(:,:,T1+1)=sigma;
           c1_serial_save(:,:,T1+1)=c1;
           c2_serial_save(:,:,T1+1)=c2;
           
           T1=T1+1
           figure(2)
           
           plot(1:length(R),[R;Y_train'],'LineWidth',1.0)
           title(['�����ģΪ',num2str(rule_num),',�����ģ��Ϊ',num2str(i),...
               '�����з�����',num2str(T1),'��ѵ��ʵ��ֵ��Ԥ��ֵͼ��','(\tau =',num2str(tao),')'])
           xlabel('t')
           ylabel('s(t)')
           legend('Ԥ��ֵ','ʵ��ֵ')          
    end    

    wz_serial=find(RMSE_serial_save==min(RMSE_serial_save));
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

%     figure(4)
%     plot(i*data_size+1+length(R):i*data_size+length(Rs_serial)+length(R),[Rs_serial;Y_test'],'LineWidth',1.0)
%     xlabel('t')
%     ylabel('s(t)')
%     title(['���з���ʵ��ֵ������Ԥ��ֵͼ��','(\tau =',num2str(tao),')'])
%     legend('Ԥ��ֵ','ʵ��ֵ')    
  
end
    SAVE_RMSE_serial_F(rule_num,:,i)    = RMSE_serial_save;     
end

%% ���г��򲿷�
for core_num=1:CORE_NUMBER     %��������ѭ��
    
    parpool(core_num);          %��MATLAB�п���������Ϊcore_num�Ĳ��м����  
    
for rule_num=1:RULR_NUMBER
    
    N=RULE_WIDTH*rule_num;                     %������Ŀ��ѭ��

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
           [R1,R2,R(i11)]=sfls_type2(X_train(i11,:),M1,M2,sigma,c1,c2); %R(i)=1/2(R1+R2)
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
    while  T2 <Tmax && RMSE >=setRMSE_accept
        
           sigma_help=zeros(N,n,core_num);
           M1_help=zeros(N,n,core_num);
           M2_help=zeros(N,n,core_num);
           c1_help=zeros(N,1,core_num);
           c2_help=zeros(N,1,core_num);  
           sigma22=sigma;M22=M2;M12=M1;
           c12=c1;c22=c2;   
           XX_train=zeros(circulation_num,n,core_num);
           YY_train=zeros(circulation_num,1,core_num); 

           tic
           parfor i22=1:core_num
                  XX_train(:,:,i22)=X_train(((i22-1)*circulation_num+1):(i22*circulation_num),:);
                  YY_train(:,:,i22)=Y_train(((i22-1)*circulation_num+1):(i22*circulation_num));
                  [M1_help(:,:,i22),M2_help(:,:,i22),c1_help(:,:,i22),c2_help(:,:,i22),sigma_help(:,:,i22)]...
                            =train_sfls_type2(XX_train(:,:,i22),YY_train(:,:,i22),M1,M2,sigma,c1,c2,alpha);
           end
           ttime21=toc
           time21=time21+ttime21
           
           for ii=2:core_num
               M1_help(:,:,1)=M1_help(:,:,1)+M1_help(:,:,ii);
               M2_help(:,:,1)=M2_help(:,:,1)+M2_help(:,:,ii);
               sigma_help(:,:,1)=sigma_help(:,:,1)+sigma_help(:,:,ii);
               c1_help(:,:,1)=c1_help(:,:,1)+c1_help(:,:,ii);
               c2_help(:,:,1)=c2_help(:,:,1)+c2_help(:,:,ii);
           end     

           sigma=sigma_help(:,:,1)/core_num;
           M1=M1_help(:,:,1)/core_num;
           M2=M2_help(:,:,1)/core_num;
           c1=c1_help(:,:,1)/core_num;
           c2=c2_help(:,:,1)/core_num;   
            
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
           time21=time21+ttime22;
     
           RMSE = sqrt(sum((R'-Y_train).^2)/length(R))
           RMSE_parallel_save(T2+1)=RMSE;
      
           T2 = T2+1
           figure(5)
                  
            plot(1:length(R),[R;Y_train'],'LineWidth',1.0)
           title(['������Ϊ',num2str(core_num),'�����ģΪ',num2str(rule_num),',�����ģ��Ϊ',num2str(i),...
               '�����з�����',num2str(T2),'��ѵ��ʵ��ֵ��Ԥ��ֵͼ��','(\tau =',num2str(tao),')'])
           xlabel('t')
           ylabel('s(t)')
           legend('Ԥ��ֵ','ʵ��ֵ')
    end
    parallel_interior_time=parallel_interior_time+time21+time22;
    
    wz=find(RMSE_parallel_save==min(RMSE_parallel_save));
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
    SAVEparallel_total_time_A(core_num,i,rule_num)    = parallel_total_time;      
    SAVEparallel_interior_time_B(core_num,i,rule_num) = parallel_interior_time;         
    SAVEspeedUP_ratio_D(core_num,i,rule_num)          = SAVEserial_total_time_C(rule_num,i)/parallel_total_time;
    SAVEparallel_efficiency_E(core_num,i,rule_num)    = SAVEserial_total_time_C(rule_num,i)/(core_num*parallel_total_time);
    SAVESserial_proportion_H(core_num,i,rule_num)     = (parallel_total_time - parallel_interior_time)/parallel_total_time;   
end       
    SAVE_RMSE__parallel_G(core_num,:,i,rule_num) = RMSE_parallel_save;
end  
   
   delete(gcp('nocreate')) %ÿһ��ѭ����ر�MATLAB���м���أ���һ�ֿ�����ͬ�������Ĳ��м����
end

%% ����������Ҫ�����ݵ�SAVE_DATA.mat�ļ��У��Ա����б�����ͼ
save('SAVEParallel_performance_parameter.mat','SAVEparallel_total_time_A','SAVEparallel_interior_time_B',...
    'SAVEserial_total_time_C','SAVEspeedUP_ratio_D','SAVEparallel_efficiency_E',...
    'SAVE_RMSE_serial_F','SAVE_RMSE__parallel_G','SAVESserial_proportion_H')
save('SAVEdata_parallel_parameter','SAVE_parallel_M1','SAVE_parallel_M2','SAVE_parallel_sigma',...
    'SAVE_parallel_c1','SAVE_parallel_c2','SAVE_parallel_Ytest','SAVE_parallel_Ypridect')
save('SAVEdata_serial_parameter','SAVE_serial_M1','SAVE_serial_M2','SAVE_serial_sigma',...
    'SAVE_serial_c1','SAVE_serial_c2','SAVE_serial_Ytest','SAVE_serial_Ypridect')



