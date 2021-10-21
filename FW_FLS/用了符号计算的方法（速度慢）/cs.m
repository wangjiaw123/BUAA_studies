clear,clc
syms x y z 
f1=cos(x)*sin(y)*tan(z)+x^2*y^3*z^4
f2=jacobian(f1,[x,y,z])
f3=char(f2)
fid=fopen('test2.txt','w');
fprintf(fid,'%s/n',f3);
fclose(fid);