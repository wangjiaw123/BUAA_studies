close all
clear,clc;
tao1=13;
tao2=18;
N=20000;%时隙为0.1，所以为0-2000s内的时间序列
[y1,t1]=Mackey_Glass(N,tao1);
[y2,t2]=Mackey_Glass(N,tao2);

figure(1)

% tao1  plot 1000-1500s之间的Mackey_Glass序列
subplot(2,2,1)
plot(t1(10000:15000),y1(10000:15000),'LineWidth',1.0);
xlabel('t')
xlabel('t','fontsize',20,'fontname','times?new?roman','FontAngle','italic');
ylabel('x(t)','fontsize',20,'fontname','times?new?roman','FontAngle','italic');
axis([10000 15001 0.58 1.24])
set(gca,'FontName','Times New Roman','FontSize',15);
grid on

% tao2  plot 1000-1500s之间的Mackey_Glass序列
subplot(2,2,2)
plot(t2(10000:15000),y2(10000:15000),'LineWidth',1.0);
xlabel('t')
xlabel('t','fontsize',20,'fontname','times?new?roman','FontAngle','italic');
ylabel('x(t)','fontsize',20,'fontname','times?new?roman','FontAngle','italic');
axis([10000 15001 0.2 1.4])
set(gca,'FontName','Times New Roman','FontSize',15);
grid on

% tao1相图时差 
subplot(2,2,3)
Dt=tao1;  
y11=y1(10000:15000);
plot(y11((Dt*10+1):end),y11(1:end-10*Dt),'LineWidth',0.5);
xlabel('s(t)','fontsize',20,'fontname','times?new?roman','FontAngle','italic');
ylabel('s(t-d)','fontsize',20,'fontname','times?new?roman','FontAngle','italic');
set(gca,'FontName','Times New Roman','FontSize',15);
grid on

% tao2相图时差 
subplot(2,2,4)
Dt=tao2;  
y22=y2(10000:15000);
plot(y22((Dt*10+1):end),y22(1:end-10*Dt),'LineWidth',0.5);
xlabel('s(t)','fontsize',20,'fontname','times?new?roman','FontAngle','italic');
ylabel('s(t-d)','fontsize',20,'fontname','times?new?roman','FontAngle','italic');
set(gca,'FontName','Times New Roman','FontSize',15);
grid on

%plot tao2 1000-2000s之间的Mackey_Glass序列
figure(2)
plot(t2(10000:20000),y2(10000:20000),'LineWidth',1.0);
set(gca,'FontName','Times New Roman','FontSize',15);
