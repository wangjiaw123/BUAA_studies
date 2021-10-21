function [y1,y2]=compute_gauss2(xx,m1,m2,sigma)
s=sigma;

if xx<m1
    y1=exp(-1/2.*((xx-m1)./(s)).^2);
elseif xx>m2
       y1=exp(-1/2.*((xx-m2)./(s)).^2); 
else
    y1=1;
end
if xx<(m1+m2)/2
    y2=exp(-1/2.*((xx-m2)./(s)).^2); 
else
   y2=exp(-1/2.*((xx-m1)./(s)).^2); 
end
end
