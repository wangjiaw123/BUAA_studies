%% 利用change_variable_method(变量变换法)求解
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
value_min_change = -value_min_change1
value_max_change
