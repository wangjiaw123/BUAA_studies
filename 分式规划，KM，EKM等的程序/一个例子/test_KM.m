
%% 利用interval_wtdavg法求解
c=(Z_up+Z_low)./2;
s=Z_up-c;
h=(W_up+W_low)./2;
delta=W_up-h;
tic
[value_min_inter,value_max_inter] = interval_wtdavg(c,s,h,delta)
toc