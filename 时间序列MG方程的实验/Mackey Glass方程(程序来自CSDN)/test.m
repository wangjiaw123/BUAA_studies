close all
clear
clc
tao=30;

Dt=tao;
N=2000;
[y,t]=Mackey_Glass(N,tao);
subplot(1,2,1)% ����ֱ�� 
plot(t,y,'LineWidth',1.0);
subplot(1,2,2)% ��ͼʱ�� 
plot(y((Dt*10+1):end),y(1:end-10*Dt),'LineWidth',0.5);

%%��Ȩ����������ΪCSDN������ddpicc0lo����ԭ�����£���ѭ CC 4.0 BY-SA ��ȨЭ�飬
%%ת���븽��ԭ�ĳ������Ӽ���������
%%ԭ�����ӣ�https://blog.csdn.net/ddpiccolo/article/details/89464435