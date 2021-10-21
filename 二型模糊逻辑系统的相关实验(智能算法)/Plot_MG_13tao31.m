clc
clear all 
close all



tao=13;
Dt=tao;
N=40020;                  %设置点数N，Mackey_Glass默认时隙间隔为0.1;
[y,t]=Mackey_Glass(N,tao);%龙格-库塔(Runge-Kutta)方法求解Mackey Glass方程
%% 原始图像
figure(1)
subplot(1,2,1)% 序列直出 
plot(t(1:10:5000),y(1:10:5000),'LineWidth',1.0);
xlabel('t')
ylabel('s(t)')
title('序列直出 ');

Y_y=y(1:10:5000);
subplot(1,2,2)% 相图时差 
plot(Y_y(tao+1:end),Y_y(1:end-tao),'LineWidth',0.5)
xlabel('s(t)')
ylabel('s(t-\tau)')
title('相图时差');



tao=31;
Dt=tao;
N=40020;                  %设置点数N，Mackey_Glass默认时隙间隔为0.1;
[y,t]=Mackey_Glass(N,tao);%龙格-库塔(Runge-Kutta)方法求解Mackey Glass方程
%% 原始图像
figure(2)
subplot(1,2,1)% 序列直出 
plot(t(1:10:5000),y(1:10:5000),'LineWidth',1.0);
xlabel('t')
ylabel('s(t)')
title('序列直出 ');

Y_y=y(1:10:10000);
subplot(1,2,2)% 相图时差 
plot(Y_y(tao+1:end),Y_y(1:end-tao),'LineWidth',0.5)
xlabel('s(t)')
ylabel('s(t-\tau)')
title('相图时差');