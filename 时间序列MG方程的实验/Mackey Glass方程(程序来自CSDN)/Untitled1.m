%%To obtain the time series value at integer points, the fourth-order Runge-Kutta 
%%method was used to find the numerical solution to the previous MG equation. 
%%It was assumed that x(0)=1.2, ��=17, and x(t)=0 for t<0.
%The result was saved in the file mgdata.dat. 
close all
clear
clc
load mgdata.dat
time = mgdata(:,1);
x = mgdata(:, 2);
figure(1)
plot(time,x)
title('Mackey-Glass Chaotic Time Series')
xlabel('Time (sec)')
ylabel('x(t)')
%��Ȩ����������ΪCSDN������ddpicc0lo����ԭ�����£���ѭ CC 4.0 BY-SA ��ȨЭ�飬ת���븽��ԭ�ĳ������Ӽ���������
%ԭ�����ӣ�https://blog.csdn.net/ddpiccolo/article/details/89464435