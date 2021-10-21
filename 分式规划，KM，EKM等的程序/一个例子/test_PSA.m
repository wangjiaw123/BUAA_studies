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
value_min_primal = -F2
value_max_primal