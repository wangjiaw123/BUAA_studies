%%To obtain the time series value at integer points, the fourth-order Runge-Kutta 
%%method was used to find the numerical solution to the previous MG equation. 
%%It was assumed that x(0)=1.2, τ=17, and x(t)=0 for t<0.
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
%版权声明：本文为CSDN博主「ddpicc0lo」的原创文章，遵循 CC 4.0 BY-SA 版权协议，转载请附上原文出处链接及本声明。
%原文链接：https://blog.csdn.net/ddpiccolo/article/details/89464435