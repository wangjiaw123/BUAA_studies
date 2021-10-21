clc
clear all 
close all



tao=13;
Dt=tao;
N=40020;                  %���õ���N��Mackey_GlassĬ��ʱ϶���Ϊ0.1;
[y,t]=Mackey_Glass(N,tao);%����-����(Runge-Kutta)�������Mackey Glass����
%% ԭʼͼ��
figure(1)
subplot(1,2,1)% ����ֱ�� 
plot(t(1:10:5000),y(1:10:5000),'LineWidth',1.0);
xlabel('t')
ylabel('s(t)')
title('����ֱ�� ');

Y_y=y(1:10:5000);
subplot(1,2,2)% ��ͼʱ�� 
plot(Y_y(tao+1:end),Y_y(1:end-tao),'LineWidth',0.5)
xlabel('s(t)')
ylabel('s(t-\tau)')
title('��ͼʱ��');



tao=31;
Dt=tao;
N=40020;                  %���õ���N��Mackey_GlassĬ��ʱ϶���Ϊ0.1;
[y,t]=Mackey_Glass(N,tao);%����-����(Runge-Kutta)�������Mackey Glass����
%% ԭʼͼ��
figure(2)
subplot(1,2,1)% ����ֱ�� 
plot(t(1:10:5000),y(1:10:5000),'LineWidth',1.0);
xlabel('t')
ylabel('s(t)')
title('����ֱ�� ');

Y_y=y(1:10:10000);
subplot(1,2,2)% ��ͼʱ�� 
plot(Y_y(tao+1:end),Y_y(1:end-tao),'LineWidth',0.5)
xlabel('s(t)')
ylabel('s(t-\tau)')
title('��ͼʱ��');