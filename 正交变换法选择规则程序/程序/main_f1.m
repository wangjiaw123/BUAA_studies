close all
clear 
clc
fls_BP1

r1=20;%设置想要留下的规则数目(可以修改)
r2=15;%设置想要留下的规则数目(可以修改)
r3=10;%设置想要留下的规则数目(可以修改)
r4=5;%设置想要留下的规则数目(可以修改)
figure
r=r1;
hold on
OLS_method
r=r2;
OLS_method
r=r3;
OLS_method
r=r4;
OLS_method

himage=findobj('tag','axes1');
set(himage,'visible','off');
set(himage,'handlevisibility','off');
hold off

figure
r=r1;
hold on
ED_method
r=r2;
ED_method
r=r3;
ED_method
r=r4;
ED_method
hold off

figure
r=r1;
hold on
SVD_QR_method
r=r2;
SVD_QR_method
r=r3;
SVD_QR_method
r=r4;
SVD_QR_method
hold off

figure
r=r1;
hold on
TLS_method
r=r2;
TLS_method
r=r3;
TLS_method
r=r4;
TLS_method
hold off

figure
r=r1;
hold on
Direct_SVDmethod
r=r2;
Direct_SVDmethod
r=r3;
Direct_SVDmethod
r=r4;
Direct_SVDmethod
hold off

% figure
% plot(value_BP1)
% hold on
% plot(value_DSVD)
% plot(value_OLS)
% plot(value_SVD_QR)
% plot(value_TLS)
% legend('value_BP1','value_DSVD','value_OLS','value_SVD_QR','value_TLS')
% 
% figure
% plot(e_BP1)
% hold on
% plot(e_DSVD)
% plot(e_OLS)
% plot(e_SVD_QR)
% plot(e_TLS)
% legend('e_BP1','e_DSVD','e_OLS','e_SVD_QR','e_TLS')
% PT=[value_BP1',value_DSVD,value_OLS,value_SVD_QR,value_TLS];
% figure
% plot(PT)
% legend('value_BP1','value_DSVD','value_OLS','value_SVD_QR','value_TLS')