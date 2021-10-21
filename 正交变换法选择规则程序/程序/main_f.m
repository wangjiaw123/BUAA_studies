close all
clear 
clc
fls_BP
r=25;%设置想要留下的规则数目(可以修改)
OLS_method
ED_method
SVD_QR_method
TLS_method
Direct_SVDmethod
figure
plot(value_BP)
hold on
plot(value_DSVD)
plot(value_OLS)
plot(value_SVD_QR)
plot(value_TLS)
legend('value_BP1','value_DSVD','value_OLS','value_SVD_QR','value_TLS')

figure
plot(e_BP)
hold on
plot(e_DSVD)
plot(e_OLS)
plot(e_SVD_QR)
plot(e_TLS)
legend('e_BP1','e_DSVD','e_OLS','e_SVD_QR','e_TLS')