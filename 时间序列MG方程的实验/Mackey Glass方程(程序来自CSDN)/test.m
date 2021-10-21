close all
clear
clc
tao=30;

Dt=tao;
N=2000;
[y,t]=Mackey_Glass(N,tao);
subplot(1,2,1)% 序列直出 
plot(t,y,'LineWidth',1.0);
subplot(1,2,2)% 相图时差 
plot(y((Dt*10+1):end),y(1:end-10*Dt),'LineWidth',0.5);

%%版权声明：本文为CSDN博主「ddpicc0lo」的原创文章，遵循 CC 4.0 BY-SA 版权协议，
%%转载请附上原文出处链接及本声明。
%%原文链接：https://blog.csdn.net/ddpiccolo/article/details/89464435