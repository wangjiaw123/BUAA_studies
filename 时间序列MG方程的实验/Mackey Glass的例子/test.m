close all
clear
clc
%% ��������
star_time=clock;
tao=31;
Dt=tao;
N=20020;%���õ���N��ʱ����Ϊ0.1;
[y,t]=Mackey_Glass(N,tao);%����-����(Runge-Kutta)�������Mackey Glass����
figure(1)
subplot(1,2,1)% ����ֱ�� 
plot(t,y,'LineWidth',1.0);
xlabel('t')
ylabel('s(t)')
title('����ֱ�� ');
subplot(1,2,2)% ��ͼʱ�� 
plot(y((Dt*10+1+10000):end),y(10001:end-10*Dt),'LineWidth',0.5);
xlabel('s(t)')
ylabel('s(t-\tau)')
title('��ͼʱ��');

n_tr = 1000;
x_star = zeros(1,n_tr);
for i = 1:n_tr
    x_star(i) = y(10001+i*10);
end
%% ѵ�������Լ��������ݶ�
X_train = [x_star(1:500);x_star(2:501);x_star(3:502);x_star(4:503);]';
Y_train = x_star(5:504)';
X_test = [x_star(505:996);x_star(506:997);x_star(507:998);x_star(508:999);]';
Y_test = x_star(509:1000)';

for i=1:(1500/500)  %�����������������벢�г������Ա�
    X_train=[X_train;X_train];
    Y_train=[Y_train;Y_train];
end
%% 
n = length(X_train(1,:));
N = 20;%���ù�����Ŀ20
M1 = rand(N,n);
M2 = M1 + 0.3;0.3
sigma = rand(N,n);
c1 = rand(N,1);
c2 = c1+0.55;0.6/0.55
cc1 = c1;
cc2 = c2;
alpha =0.4;%KMʱȡ0.4
Err = 0.2;
[R1,R2,R]=sfls_type2(X_train,M1,M2,sigma,c1,c2);
Error = sqrt(sum((R'-Y_train).^2))
T = 1;
Tmax = 20;

M1_save=zeros(N,n,Tmax);
M2_save=zeros(N,n,Tmax);
c1_save=zeros(N,1,Tmax);
c2_save=zeros(N,1,Tmax);
sigma_save=zeros(N,n,Tmax);
Error_save=[];
R_save=zeros(Tmax,length(Y_train));
total_train_time=0;
while  T <=Tmax %����ѵ������
       %alpha=1/T
       Error1 = Error;
      %[M1,M2,c1,c2,sigma]=train_sfls(X_train(:,:),Y_train(:,:),M1,M2,sigma,c1,c2,alpha);
      %[M1,M2,c1,c2,sigma,I2l,I2u,I1u,I1l]=train_sfls_type2(X_train,Y_train,M1,M2,sigma,c1,c2,alpha);
      tic
      [M1,M2,c1,c2,sigma,I2l,I2u,I1u,I1l]=train_sfls_type2(X_train,Y_train,M1,M2,sigma,c1,c2,alpha);
      p_time=toc
      total_train_time=total_train_time+p_time;
      M1_save(:,:,T)=M1;
      M2_save(:,:,T)=M2;
      c1_save(:,:,T)=c1;
      c2_save(:,:,T)=c2;
      sigma_save(:,:,T)=sigma;

      [R1,R2,R]=sfls_type2(X_train,M1,M2,sigma,c1,c2);  
      R_save(T,:)=R;
      RMSE = sqrt(sum((R'-Y_train).^2)/length(R))
      Error = sqrt(sum((R'-Y_train).^2))
      Error_save=[Error_save Error];
      T = T+1
      figure(2)
      plot([1000+1:1000+length(R)],[R;Y_train'],'LineWidth',1.0)
      title(['��',num2str(T),'��ѵ��ʵ��ֵ��Ԥ��ֵͼ��','(\tau =',num2str(tao),')'])
      xlabel('t')
      ylabel('s(t)')
      legend('Ԥ��ֵ','ʵ��ֵ')
end
wz=find(Error_save==min(Error_save));
[R11,R22,Rs]=sfls_type2(X_test,M1_save(:,:,wz),M2_save(:,:,wz),sigma_save(:,:,wz),c1_save(:,:,wz),c2_save(:,:,wz));
RMSE = sqrt(sum((Rs'-Y_test).^2)/length(Rs))
figure(5)
plot([1000+1:1000+length(R_save(1,:))],[R_save(wz,:);Y_train'],'LineWidth',1.0)
title(['��',num2str(wz),'��ѵ��ʵ��ֵ��Ԥ��ֵͼ��','(\tau =',num2str(tao),')'])
xlabel('t')
ylabel('s(t)')
legend('Ԥ��ֵ','ʵ��ֵ')

figure(3)
subplot(1,2,1)
plot([1000+1:1000+length(Rs)],Rs,'LineWidth',1.0)
xlabel('t')
ylabel('s(t)')
title('����ֱ�� ');
subplot(1,2,2)
plot(Rs((Dt*10+1):end),Rs(1:end-10*Dt),'LineWidth',0.5)
xlabel('s(t)')
ylabel('s(t-\tau)')
title('��ͼʱ��');
figure(4)
plot([1001+length(R):1000+length(Rs)+length(R)],[Rs;Y_test'],'LineWidth',1.0)
xlabel('t')
ylabel('s(t)')
title(['ʵ��ֵ������Ԥ��ֵͼ��','(\tau =',num2str(tao),')'])
legend('Ԥ��ֵ','ʵ��ֵ')
total_train_time
end_time=clock;
serial_total_time=etime(end_time,star_time)
















