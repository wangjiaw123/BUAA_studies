close all
clear 
clc
fls_BP1
hold on
r=25;%������Ҫ���µĹ�����Ŀ(�����޸�)
figure
OLS_method
figure
ED_method
figure
SVD_QR_method
figure
TLS_method
figure
Direct_SVDmethod

figure
hold on
plot(value_DSVD)
plot(value_OLS)
plot(value_SVD_QR)
plot(value_TLS)
legend('value_BP1','value_DSVD','value_OLS','value_SVD_QR','value_TLS')

figure
plot(e_BP1)
hold on
plot(e_DSVD)
plot(e_OLS)
plot(e_SVD_QR)
plot(e_TLS)
legend('e_BP1','e_DSVD','e_OLS','e_SVD_QR','e_TLS')