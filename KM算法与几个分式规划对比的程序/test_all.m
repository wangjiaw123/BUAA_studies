close all
clear
clc
%% 生成测试数据
nn=2500;
data_gap=2;
data_num=nn/data_gap;
ZZ_LOW = 10*rand(nn,1);
ZZ_UP = ZZ_LOW+10*rand(nn,1);
WW_LOW = rand(nn,1);
WW_UP = WW_LOW+rand(nn,1);
%%
MonteCarlo=1;
TIME=zeros(7,data_num,MonteCarlo);
TIME_MEAN=zeros(7,data_num);

%% 
for MC =1:MonteCarlo  %蒙特卡罗,似乎不需要
    for dt=1:data_num
        dt*data_gap
        Z_low=ZZ_LOW(1:dt*data_gap);
        Z_up=ZZ_UP(1:dt*data_gap);
        W_low=WW_LOW(1:dt*data_gap);
        W_up=WW_UP(1:dt*data_gap);
        
        %% EKM
        disp('EKM')
        tic
        [EKM_min_x,EKM_min,EKM_max_x,EKM_max] = EKM(Z_low,Z_up,W_low,W_up);
        EKM_time=toc;
        TIME(1,dt,MC)=EKM_time;
        %% DA
%         disp('DA')        
%         tic
%         [DA_min_x,DA_min,DA_max_x,DA_max] = DA(Z_low,Z_up,W_low,W_up);
%         DA_time=toc;
%         TIME(2,dt,MC)=DA_time;        
        %% KM
        disp('KM')        
        c=(Z_up+Z_low)./2;
        s=Z_up-c;
        h=(W_up+W_low)./2;
        delta=W_up-h;
        tic
        [KM_min,KM_max] = interval_wtdavg(c,s,h,delta);
        KM_time=toc;   
        TIME(3,dt,MC)=KM_time;
%        %% 利用iterative_method(迭代算法)
%         disp('iterative_method')        
%         n_psa = length(W_up);
%         A = zeros(2*n_psa,3*n_psa);
%         A(1:n_psa,1:n_psa) = eye(n_psa);
%         A(n_psa+1:2*n_psa,1:n_psa) = -eye(n_psa);
%         A(1:2*n_psa,n_psa+1:3*n_psa) = eye(2*n_psa);
%         b = [W_up;-W_low];
%         d1 = [ones(1,n_psa) zeros(1,2*n_psa)];
%         d0 = 0;
%         % d1 =zeros(1,3*n_psa);
%         % d0 = 1;
%         c0 = 0;
%         c1 = [Z_up' zeros(1,2*n_psa)]; %设置c1求最大值
%         c2 = [Z_low' zeros(1,2*n_psa)];%设置c2求最小值
%         tic
%         [ X_max_iterative,value_max_iterative ] = iterative_method( c1,c0,d1,d0,A,b );
%         [ X_min_iterative,value_min_iterative1 ] = iterative_method( -c2,-c0,d1,d0,A,b );
%         iterative_time=toc;
%         value_min_iterative = -value_min_iterative1;
%         value_max_iterative;
%         TIME(4,dt,MC)=iterative_time;           
%         %% Primal simplex algorithm(主单纯型算法) 
%         disp('Primal simplex algorithm')        
%         n_psa = length(W_up);
%         A = zeros(2*n_psa,3*n_psa);
%         A(1:n_psa,1:n_psa) = eye(n_psa);
%         A(n_psa+1:2*n_psa,1:n_psa) = -eye(n_psa);
%         A(1:2*n_psa,n_psa+1:3*n_psa) = eye(2*n_psa);
%         b = [W_up;-W_low];
%         d1 = [ones(1,n_psa) zeros(1,2*n_psa)];
%         d0 = 0;
%         % d1 =zeros(1,3*n_psa);
%         % d0 = 1;
%         c0 = 0;
%         c1 = [Z_up' zeros(1,2*n_psa)]; %设置c1求最大值
%         c2 = [Z_low' zeros(1,2*n_psa)];%设置c2求最小值
%         tic
%         [ sheet1,sheet2,cow_location1,value_max_primal,cow_location_save1,T1] = primal_simplex_algorithm( c1,c0,d1,d0,A,b );
%         [ sheet11,sheet22,cow_location2,F2,cow_location_save2,T2 ] = primal_simplex_algorithm( -c2,-c0,d1,d0,A,b );
%         PSA_time=toc;
%         PSA_value_min_primal = -F2     
%         TIME(5,dt,MC)=PSA_time;
%         %% approximate_method(邻近算法)
%         disp('approximate_method')      
%         n_psa = length(W_up);
%         A = zeros(2*n_psa,3*n_psa);
%         A(1:n_psa,1:n_psa) = eye(n_psa);
%         A(n_psa+1:2*n_psa,1:n_psa) = -eye(n_psa);
%         A(1:2*n_psa,n_psa+1:3*n_psa) = eye(2*n_psa);
%         b = [W_up;-W_low];
%         d1 = [ones(1,n_psa) zeros(1,2*n_psa)];
%         d0 = 0;
%         % d1 =zeros(1,3*n_psa);
%         % d0 = 1;
%         c0 = 0;
%         c1 = [Z_up' zeros(1,2*n_psa)]; %设置c1求最大值
%         c2 = [Z_low' zeros(1,2*n_psa)];%设置c2求最小值
%         tic
%         [ X_max_approximate,value_max_approximate ] = approximate_method( c1,d1,A,b );
%         [ X_min_approximate,value_min_approximate1 ] = approximate_method( -c2,d1,A,b );
%         approximate_time=toc;
%         value_min_approximate = -value_min_approximate1;
%         value_max_approximate;
%         TIME(6,dt,MC)=approximate_time;        
%       
%         %% 利用change_variable_method(变量变换法)求解
%         disp('change_variable_method')          
%         n_psa = length(W_up);
%         A = zeros(2*n_psa,3*n_psa);
%         A(1:n_psa,1:n_psa) = eye(n_psa);
%         A(n_psa+1:2*n_psa,1:n_psa) = -eye(n_psa);
%         A(1:2*n_psa,n_psa+1:3*n_psa) = eye(2*n_psa);
%         b = [W_up;-W_low];
%         d1 = [ones(1,n_psa) zeros(1,2*n_psa)];
%         d0 = 0;
%         % d1 =zeros(1,3*n_psa);
%         % d0 = 1;
%         c0 = 0;
%         c1 = [Z_up' zeros(1,2*n_psa)]; %设置c1求最大值
%         c2 = [Z_low' zeros(1,2*n_psa)];%设置c2求最小值
%         AA = zeros(2*n_psa,n_psa);
%         AA(1:n_psa,1:n_psa) = eye(n_psa);
%         AA(n_psa+1:2*n_psa,1:n_psa) = -eye(n_psa);
%         c11 = Z_up';
%         d11 = ones(1,n_psa);
%         % d11 = zeros(1,n_psa);
%         tic
%         [ X_max_change,value_max_change ] = change_variable_method( c11,c0,d11,d0,AA,b );
%         c22 = Z_low';
%         [ X_min_change,value_min_change1 ] = change_variable_method( -c22,-c0,d11,d0,AA,b );
%         change_variable_time=toc;
%         value_min_change = -value_min_change1;
%         value_max_change;
%         TIME(7,dt,MC)=change_variable_time;         
        
    end    
end

for MCMC=1:MonteCarlo
    TIME_MEAN=TIME_MEAN+TIME(:,:,MCMC);
end
TIME_MEAN=TIME_MEAN/MonteCarlo;

RGB=[0,0,0; %黑色
    1,0,0;  %红色
    1,0.5,0;
    1,0.5,0.5;
    1,0,0.5;
    0,1,0.5;
    1,0,1;
    1,1,0;
    0.5,1,0;
    0,0,1;
    0.5,0.2,1;
    1,0,1;
    0.5,0.5,0.5;];


set(0,'defaultfigurecolor','w')

% figure(1)
% hold on
% box on
% axis([20 nn 0 0.035])
% for k=1:3
%     plot([data_gap:data_gap:nn],TIME_MEAN(k,:),'Color',RGB(k,:))
% end
% plot([data_gap:data_gap:nn],TIME_MEAN(4,:),'Color',RGB(4,:))
% title(['EKM,DA,KM,IM算法的计算时间比较'])
% ylabel('time')
% xlabel('n')
% legend('EKM','DA','KM','IM','Location', 'northeastoutside' )
% hold off

% figure(2)
% hold on
% box on
% axis([20 nn 0 1.25])
% for k=4:7
%     plot([data_gap:data_gap:nn],TIME_MEAN(k,:),'Color',RGB(k,:))
% end
% title(['IM,PSA,AM,CVM算法的计算时间比较'])
% ylabel('time')
% xlabel('n')
% legend('IM','PSA','AM','CVM','Location', 'northeastoutside' )
% hold off

figure(3)
hold on
box on
axis([20 nn 0 0.002])
for k=1:2:3
    plot([data_gap:data_gap:nn],TIME_MEAN(k,:),'Color',RGB(k,:))
end
title(['EKM,KM算法的计算时间'])
ylabel('time')
xlabel('n')
legend('EKM','KM','Location', 'northeastoutside' )
hold off


save('DATA_SAVE','TIME_MEAN','TIME')