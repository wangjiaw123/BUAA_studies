close all
clear
clc
%% 生成数据
star_time=clock;
tao=31;
Dt=tao;
N=20020;%设置点数N，时间间隔为0.1;
[y,t]=Mackey_Glass(N,tao);%龙格-库塔(Runge-Kutta)方法求解Mackey Glass方程
figure(1)
subplot(1,2,1)% 序列直出 
plot(t,y,'LineWidth',1.0);
xlabel('t')
ylabel('s(t)')
title('序列直出 ');
subplot(1,2,2)% 相图时差 
plot(y((Dt*10+1+10000):end),y(10001:end-10*Dt),'LineWidth',0.5);
xlabel('s(t)')
ylabel('s(t-\tau)')
title('相图时差');

n_tr = 1000;
x_star = zeros(1,n_tr);
for i = 1:n_tr
    x_star(i) = y(10001+i*10);
end
%% 训练数据以及测试数据对
X_train = [x_star(1:500);x_star(2:501);x_star(3:502);x_star(4:503);]';
Y_train = x_star(5:504)';

% XXX_train=X_train;
% YYY_train=Y_train;
% for i=1:1 %用于增大数据量，与并行程序作对比
%     XXX_train=[XXX_train;X_train];
%     YYY_train=[YYY_train;Y_train];
% end
% X_train=XXX_train;
% Y_train=YYY_train;

X_test = [x_star(505:996);x_star(506:997);x_star(507:998);x_star(508:999);]';
Y_test = x_star(509:1000)';
%% 
n = length(X_train(1,:));
N = 20;%设置规则数目20
M1 = rand(N,n);
M2 = M1 + 0.3;
sigma = rand(N,n);
c1 = rand(N,1);
c2 = c1+0.55;
cc1 = c1;
cc2 = c2;
alpha =0.2;%KM时取0.4
Err = 0.2;
R=zeros(1,n);
tic
parfor i=1:length(X_train(:,1))
[R1,R2,R(i)]=sfls_type2(X_train(i,:),M1,M2,sigma,c1,c2);
end
compute_time1=toc
Error = sqrt(sum((R'-Y_train).^2));
T = 1;
Tmax = 50;

M1_save=zeros(N,n,Tmax);
M2_save=zeros(N,n,Tmax);
c1_save=zeros(N,1,Tmax);
c2_save=zeros(N,1,Tmax);
sigma_save=zeros(N,n,Tmax);
Error_save=[];


core_num=6;%设置核心数
row_num=length(X_train(:,1));
%parpool(core_num)
lack_num=core_num-mod(row_num,core_num);
R_save=zeros(Tmax,length(Y_train)+lack_num);
for k=1:lack_num
    X_train=[X_train' X_train(k,:)']';
    Y_train=[Y_train' Y_train(k,:)']';
end
total_num=row_num+lack_num;
circulation_num=total_num/core_num;
XX_train=zeros(core_num,length(X_train));
paraller_train_time=0;
compute_time2=0
while  T <=Tmax %设置训练轮数
       Error1 = Error;
      %[M1,M2,c1,c2,sigma]=train_sfls(X_train(:,:),Y_train(:,:),M1,M2,sigma,c1,c2,alpha);
      %[M1,M2,c1,c2,sigma,I2l,I2u,I1u,I1l]=train_sfls_type2(X_train,Y_train,M1,M2,sigma,c1,c2,alpha);
      
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
          parfor i=1:core_num
          XX_train(:,:,i)=X_train(((i-1)*circulation_num+1):(i*circulation_num),:);
          YY_train(:,:,i)=Y_train(((i-1)*circulation_num+1):(i*circulation_num));
          [M1_help(:,:,i),M2_help(:,:,i),c1_help(:,:,i),c2_help(:,:,i),sigma_help(:,:,i)]...
              =train_sfls_type2(XX_train(:,:,i),YY_train(:,:,i),M1,M2,sigma,c1,c2,alpha);
          end
          for ii=2:core_num
              M1_help(:,:,1)=M1_help(:,:,1)+M1_help(:,:,ii);
              M2_help(:,:,1)=M2_help(:,:,1)+M2_help(:,:,ii);
              sigma_help(:,:,1)=sigma_help(:,:,1)+sigma_help(:,:,ii);
              c1_help(:,:,1)=c1_help(:,:,1)+c1_help(:,:,ii);
              c2_help(:,:,1)=c2_help(:,:,1)+c2_help(:,:,ii);
          end     
          sigma=sigma22+(sigma_help(:,:,1)-core_num*sigma22)/core_num;
          M1=M12+(M1_help(:,:,1)-core_num*M12)/core_num;
          M2=M22+(M2_help(:,:,1)-core_num*M22)/core_num;
          c1=c12+(c1_help(:,:,1)-core_num*c12)/core_num;
          c2=c22+(c2_help(:,:,1)-core_num*c22)/core_num;      
          par_time=toc
          paraller_train_time=paraller_train_time+par_time;          

%       tic
%       for k=1:circulation_num-1
%           XX_train=X_train((k-1)*core_num+1:k*core_num,:);
%           YY_train=Y_train((k-1)*core_num+1:k*core_num);
%           sigma_help=zeros(N,n);
%           M1_help=zeros(N,n);
%           M2_help=zeros(N,n);
%           c1_help=zeros(N,1);
%           c2_help=zeros(N,1);
%           sigma22=sigma;M22=M2;M12=M1;
%           c12=c1;c22=c2;
%           spmd
%           [M1,M2,c1,c2,sigma]=train_sfls_type2(XX_train(labindex,:),YY_train(labindex),M1,M2,sigma,c1,c2,alpha);
%           end
% %           for ii=1:core_num
% %               M1_help=M1_help+cell2mat(M1(ii));
% %               M2_help=M2_help+cell2mat(M2(ii));
% %               sigma_help=sigma_help+cell2mat(sigma(ii));
% %               c1_help=c1_help+cell2mat(c1(ii));
% %               c2_help=c2_help+cell2mat(c2(ii));
% %           end
%           for ii=1:core_num
%               M1_help=M1_help+M1{ii};
%               M2_help=M2_help+M2{ii};
%               sigma_help=sigma_help+sigma{ii};
%               c1_help=c1_help+c1{ii};
%               c2_help=c2_help+c2{ii};
%           end
%           sigma=sigma22+(sigma_help-core_num*sigma22)/core_num;
%           M1=M12+(M1_help-core_num*M12)/core_num;
%           M2=M22+(M2_help-core_num*M22)/core_num;
%           c1=c12+(c1_help-core_num*c12)/core_num;
%           c2=c22+(c2_help-core_num*c22)/core_num;
%       end
%       par_time=toc
%       paraller_time=paraller_time+par_time;
      
      
%       [M1,M2,c1,c2,sigma]=train_sfls_type2(X_train,Y_train,M1,M2,sigma,c1,c2,alpha);
      M1_save(:,:,T)=M1;
      M2_save(:,:,T)=M2;
      c1_save(:,:,T)=c1;
      c2_save(:,:,T)=c2;
      sigma_save(:,:,T)=sigma;
      
      tic
      parfor i=1:length(X_train(:,1))
            [R1,R2,R(i)]=sfls_type2(X_train(i,:),M1,M2,sigma,c1,c2);
      end
      comp=toc
      compute_time2=compute_time2+comp;
      %[R1,R2,R]=sfls_type2(X_train,M1,M2,sigma,c1,c2);  
      R_save(T,:)=R;
      
      RMSE = sqrt(sum((R'-Y_train).^2)/length(R))
      Error = sqrt(sum((R'-Y_train).^2))
      Error_save=[Error_save Error];
      
      T = T+1
      figure(2)
      plot(1000+1:1000+length(R),[R;Y_train'],'LineWidth',1.0)
      title(['第',num2str(T),'次训练实际值与预测值图像','(\tau =',num2str(tao),')'])
      xlabel('t')
      ylabel('s(t)')
      legend('预测值','实际值')
end
figure(6)
plot(Error_save)
title('误差随训练次数的变化');
xlabel('k');
ylabel('train error')

wz=find(Error_save==min(Error_save));
Rs=zeros(1,length(X_test(:,1)));
tic
parfor i=1:length(X_test(:,1))
[R11,R22,Rs(i)]=sfls_type2(X_test(i,:),M1_save(:,:,wz),M2_save(:,:,wz),sigma_save(:,:,wz),c1_save(:,:,wz),c2_save(:,:,wz));
end
compute_time3=toc
RMSE = sqrt(sum((Rs'-Y_test).^2)/length(Rs))
figure(5)
plot(1000+1:1000+length(R_save(1,:)),[R_save(wz,:);Y_train'],'LineWidth',1.0)
title(['第',num2str(wz),'次训练实际值与预测值图像','(\tau =',num2str(tao),')'])
xlabel('t')
ylabel('s(t)')
legend('预测值','实际值')

figure(3)
subplot(1,2,1)
plot(1000+1:1000+length(Rs),Rs,'LineWidth',1.0)
xlabel('t')
ylabel('s(t)')
title('序列直出 ');
subplot(1,2,2)
plot(Rs((Dt*10+1):end),Rs(1:end-10*Dt),'LineWidth',0.5)
xlabel('s(t)')
ylabel('s(t-\tau)')
title('相图时差');
figure(4)
plot(1001+length(R):1000+length(Rs)+length(R),[Rs;Y_test'],'LineWidth',1.0)
xlabel('t')
ylabel('s(t)')
title(['实际值与最终预测值图像','(\tau =',num2str(tao),')'])
legend('预测值','实际值')
paraller_total_time=compute_time1+compute_time2+compute_time3+paraller_train_time
paraller_train_time
end_time=clock;
total_time=etime(end_time,star_time)














