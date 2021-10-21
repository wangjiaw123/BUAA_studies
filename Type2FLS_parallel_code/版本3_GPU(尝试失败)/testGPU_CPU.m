close all
clear
clc

% 生成数据
%star_parallel_time=clock;
tao=31;
Dt=tao;
data_size = gpuArray(200);
data_size1 = 200;
N=40020;%设置点数N，Mackey_Glass默认时隙间隔为0.1;
[y,t]=Mackey_Glass(N,tao);%龙格-库塔(Runge-Kutta)方法求解Mackey Glass方程
n=4;
M=30;
n_train = 4000;
x_star = zeros(1,n_train);
for i = 1:n_train
    x_star(i) = y(1+i*10);
end
x_star_GPU = gpuArray(x_star);
i=gpuArray(1);
X_train_GPU = [x_star_GPU(1:i*data_size-4);x_star_GPU(2:i*data_size-3);x_star_GPU(3:i*data_size-2);x_star_GPU(4:i*data_size-1);]';
Y_train_GPU = x_star(5:i*data_size)';
X_test_GPU = [x_star_GPU(i*data_size-3:(i+1)*data_size-4);x_star_GPU(i*data_size-2:(i+1)*data_size-3);...
         x_star_GPU(i*data_size-1:(i+1)*data_size-2);x_star_GPU(i*data_size:(i+1)*data_size-1)]';
Y_test_GPU = x_star_GPU(i*data_size+1:(i+1)*data_size)';

i1=1    
X_train = [x_star(1:i1*data_size1-4);x_star(2:i1*data_size1-3);x_star(3:i1*data_size1-2);x_star(4:i1*data_size1-1);]';
Y_train = x_star(5:i1*data_size1)';
X_test = [x_star(i1*data_size1-3:(i1+1)*data_size1-4);x_star(i1*data_size1-2:(i1+1)*data_size1-3);...
            x_star(i1*data_size1-1:(i1+1)*data_size1-2);x_star(i1*data_size1:(i1+1)*data_size1-1)]';
Y_test = x_star(i1*data_size1+1:(i1+1)*data_size1)';    

    N1=gather(M);
    n1=gather(n);
    M1 = gpuArray(rand(N1,n1));
    M2 = M1 + 0.8;
    sigma = gpuArray(rand(N1,n1));
    c1 = gpuArray(rand(N1,1));
    c2 = c1+1;   
 tic   
 R=sfls_type2(X_train_GPU,M1,M2,sigma,c1,c2);   
 toc 

    M11 = rand(M,n);
    M21 = M11 + 0.8;
    sigma1 = rand(N,n);
    c11 = rand(M,1);
    c21 = c11+1; 
 tic
 R1=sfls_type2(X_train,M11,M21,sigma1,c11,c21);
 toc
%  gather(R)-R1;
%  tic
%  A=rand(10000,10000);
%  b=A*A';
%  toc
%  
%  tic
%  A1=gpuArray(rand(10000,10000));
%  B1=A1*A1';
%  toc
 
 
 
 
 
 