%% 利用approximate_method(邻近算法)求解
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
[ X_max_approximate,value_max_approximate ] = approximate_method( c1,d1,A,b );
[ X_min_approximate,value_min_approximate1 ] = approximate_method( -c2,d1,A,b );
toc
value_min_approximate = -value_min_approximate1
value_max_approximate