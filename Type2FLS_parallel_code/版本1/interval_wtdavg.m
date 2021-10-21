function[l_out,r_out] = interval_wtdavg(c,s,h,delta)

lower = h - delta ;
upper = h + delta ;

l_out = adapt(c-s,lower,upper,-1) ; 
r_out = adapt(c+s,lower,upper,1) ; 

return ;

