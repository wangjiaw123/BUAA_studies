%main.m
%分别使用《基于不确定规则的模糊逻辑系统》中的程序interval_wtdavg.m以及Primal 
%    simplex algorithm(主单纯型算法),approximate_method(邻近算法),iterative
%   _method(迭代算法),change_variable_method(变量变换法)，计算sum(Z*W)/sum(W)
%   的最大值value_max和最小值value_min.
%   Z_low(i) <= Z(i) <= Z_up(i);W_low(i) <= W(i) <= W_up(i);
%   Z_low,Z_up,W_low,W_up都是列向量
close all
clear
clc
%% 一个例子的初始化（设置参数范围）
% Z_low = [-1.56,;1.3;1.5;-0.3;-0.14;5.3;4.5;2.3;5.66;8.99];
% Z_up = [3.95;6.22;2.55;1.5;3.45;9.78;8.65;9.56;13.56;9.21];
% W_low = [1.22;0.98;0.56;5.66;4.85;2.56;9.66;1.2;3.55;8.66];
% W_up = [7.85;9.65;5.68;14.56;12.24;8.34;18.65;6.33;4.62;8.70];

% Z_low = [1.5;-0.3;-2.6;0.5;3.66];
% Z_up = [2.55;1.5;2.44;1.3;5.33];
% W_low = [2.56;1.66;2.33;1.22;2.69];
% W_up = [3.68;4.56;4.56;3.21;10.28];

% Z_low = [-6.3;-1.56;1.3;-1.5];
% Z_up = [3.95;6.22;2.55;2.55];
% W_low = [1.22;0.98;4.56;3.55];
% W_up = [7.85;3.65;5.68;7.66];

nn=10;
Z_low = 100*rand(nn,1);
Z_up = Z_low+10*rand(nn,1);
W_low = 5*rand(nn,1);
W_up = Z_low+15*rand(nn,1);
%% 判断是否符合输入要求
if (length(Z_low)~=length(Z_up))||(length(Z_low)~=length(W_low))||(length(Z_low)~=length(W_up)) 
    disp('输入的列向量Z_low,Z_up,W_low,W_up长度不相同！');
end
%% 利用interval_wtdavg法求解
c=(Z_up+Z_low)./2;
s=Z_up-c;
h=(W_up+W_low)./2;
delta=W_up-h;
tic
[value_min_inter,value_max_inter] = interval_wtdavg(c,s,h,delta);
toc
%% %% 利用自己写的定理9.1的程序求解
[ value_min,X_min,value_max,X_max ] = theory9_1( Z_low,Z_up,W_low,W_up );
%% 利用Primal simplex algorithm(主单纯型算法)求解
n_psa = length(W_up);
A = zeros(2*n_psa,3*n_psa);
A(1:n_psa,1:n_psa) = eye(n_psa);
A(n_psa+1:2*n_psa,1:n_psa) = -eye(n_psa);
A(1:2*n_psa,n_psa+1:3*n_psa) = eye(2*n_psa);
b = [W_up;-W_low];
d1 = [ones(1,n_psa) zeros(1,2*n_psa)];
d0 = 0;
% d1 =zeros(1,3*n_psa);
% d0 = 1;
c0 = 0;
c1 = [Z_up' zeros(1,2*n_psa)]; %设置c1求最大值
c2 = [Z_low' zeros(1,2*n_psa)];%设置c2求最小值
tic
[ sheet1,sheet2,cow_location1,value_max_primal,cow_location_save1,T1] = primal_simplex_algorithm( c1,c0,d1,d0,A,b );
[ sheet11,sheet22,cow_location2,F2,cow_location_save2,T2 ] = primal_simplex_algorithm( -c2,-c0,d1,d0,A,b );
toc
value_min_primal = -F2;

%% 利用approximate_method(邻近算法)求解
tic
[ X_max_approximate,value_max_approximate ] = approximate_method( c1,d1,A,b );
[ X_min_approximate,value_min_approximate1 ] = approximate_method( -c2,d1,A,b );
toc
value_min_approximate = -value_min_approximate1;
%% 利用iterative_method(迭代算法)求解
tic
[ X_max_iterative,value_max_iterative ] = iterative_method( c1,c0,d1,d0,A,b );
[ X_min_iterative,value_min_iterative1 ] = iterative_method( -c2,-c0,d1,d0,A,b );
toc
value_min_iterative = -value_min_iterative1;
%% 利用change_variable_method(变量变换法)求解
AA = zeros(2*n_psa,n_psa);
AA(1:n_psa,1:n_psa) = eye(n_psa);
AA(n_psa+1:2*n_psa,1:n_psa) = -eye(n_psa);
c11 = Z_up';
d11 = ones(1,n_psa);
% d11 = zeros(1,n_psa);
tic
[ X_max_change,value_max_change ] = change_variable_method( c11,c0,d11,d0,AA,b );
c22 = Z_low';
[ X_min_change,value_min_change1 ] = change_variable_method( -c22,-c0,d11,d0,AA,b );
toc
value_min_change = -value_min_change1;
%% 输出结果
fprintf('利用课本自带的定理9.1的程序算出的最小值为：%f \n',value_min_inter)
fprintf('利用自己写的定理9.1的程序求得的最小值为：%f \n',value_min)
fprintf('利用Primal simplex algorithm(主单纯型算法)算出的最小值为：%f \n',value_min_primal)
fprintf('利用approximate_method(邻近算法)算出的最小值为：%f \n',value_min_approximate)
fprintf('利用iterative_method(迭代算法)算出的最小值为：%f \n',value_min_iterative)
fprintf('利用change_variable_method(变量变换法)算出的最小值为：%f \n',value_min_change)
fprintf('\n')
fprintf('利用课本自带的定理9.1的程序算出的最大值为：%f \n',value_max_inter)
fprintf('利用自己写的定理9.1的程序求得的最大值为：%f \n',value_max)
fprintf('利用Primal simplex algorithm(主单纯型算法)算出的最大值为：%f \n',value_max_primal)
fprintf('利用approximate_method(邻近算法)算出的最大值为：%f \n',value_max_approximate)
fprintf('利用iterative_method(迭代算法)算出的最大值为：%f \n',value_max_iterative)
fprintf('利用change_variable_method(变量变换法)算出的最大值为：%f \n',value_max_change)

% [x111,fv111,exitflag1,outpt1,lama1] = linprog(-c1,A,b,[],[],zeros(length(c1),1));
% fv111=-fv111
% [x222,fv222,exitflag1,outpt1,lama1] = linprog(c2,A,b,[],[],zeros(length(c1),1));
% fv222
