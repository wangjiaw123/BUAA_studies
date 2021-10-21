%main.m
%�ֱ�ʹ�á����ڲ�ȷ�������ģ���߼�ϵͳ���еĳ���interval_wtdavg.m�Լ�Primal 
%    simplex algorithm(���������㷨),approximate_method(�ڽ��㷨),iterative
%   _method(�����㷨),change_variable_method(�����任��)������sum(Z*W)/sum(W)
%   �����ֵvalue_max����Сֵvalue_min.
%   Z_low(i) <= Z(i) <= Z_up(i);W_low(i) <= W(i) <= W_up(i);
%   Z_low,Z_up,W_low,W_up����������
close all
clear
clc
%% һ�����ӵĳ�ʼ�������ò�����Χ��
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
%% �ж��Ƿ��������Ҫ��
if (length(Z_low)~=length(Z_up))||(length(Z_low)~=length(W_low))||(length(Z_low)~=length(W_up)) 
    disp('�����������Z_low,Z_up,W_low,W_up���Ȳ���ͬ��');
end
%% ����interval_wtdavg�����
c=(Z_up+Z_low)./2;
s=Z_up-c;
h=(W_up+W_low)./2;
delta=W_up-h;
tic
[value_min_inter,value_max_inter] = interval_wtdavg(c,s,h,delta);
toc
%% %% �����Լ�д�Ķ���9.1�ĳ������
[ value_min,X_min,value_max,X_max ] = theory9_1( Z_low,Z_up,W_low,W_up );
%% ����Primal simplex algorithm(���������㷨)���
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
c1 = [Z_up' zeros(1,2*n_psa)]; %����c1�����ֵ
c2 = [Z_low' zeros(1,2*n_psa)];%����c2����Сֵ
tic
[ sheet1,sheet2,cow_location1,value_max_primal,cow_location_save1,T1] = primal_simplex_algorithm( c1,c0,d1,d0,A,b );
[ sheet11,sheet22,cow_location2,F2,cow_location_save2,T2 ] = primal_simplex_algorithm( -c2,-c0,d1,d0,A,b );
toc
value_min_primal = -F2;

%% ����approximate_method(�ڽ��㷨)���
tic
[ X_max_approximate,value_max_approximate ] = approximate_method( c1,d1,A,b );
[ X_min_approximate,value_min_approximate1 ] = approximate_method( -c2,d1,A,b );
toc
value_min_approximate = -value_min_approximate1;
%% ����iterative_method(�����㷨)���
tic
[ X_max_iterative,value_max_iterative ] = iterative_method( c1,c0,d1,d0,A,b );
[ X_min_iterative,value_min_iterative1 ] = iterative_method( -c2,-c0,d1,d0,A,b );
toc
value_min_iterative = -value_min_iterative1;
%% ����change_variable_method(�����任��)���
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
%% ������
fprintf('���ÿα��Դ��Ķ���9.1�ĳ����������СֵΪ��%f \n',value_min_inter)
fprintf('�����Լ�д�Ķ���9.1�ĳ�����õ���СֵΪ��%f \n',value_min)
fprintf('����Primal simplex algorithm(���������㷨)�������СֵΪ��%f \n',value_min_primal)
fprintf('����approximate_method(�ڽ��㷨)�������СֵΪ��%f \n',value_min_approximate)
fprintf('����iterative_method(�����㷨)�������СֵΪ��%f \n',value_min_iterative)
fprintf('����change_variable_method(�����任��)�������СֵΪ��%f \n',value_min_change)
fprintf('\n')
fprintf('���ÿα��Դ��Ķ���9.1�ĳ�����������ֵΪ��%f \n',value_max_inter)
fprintf('�����Լ�д�Ķ���9.1�ĳ�����õ����ֵΪ��%f \n',value_max)
fprintf('����Primal simplex algorithm(���������㷨)��������ֵΪ��%f \n',value_max_primal)
fprintf('����approximate_method(�ڽ��㷨)��������ֵΪ��%f \n',value_max_approximate)
fprintf('����iterative_method(�����㷨)��������ֵΪ��%f \n',value_max_iterative)
fprintf('����change_variable_method(�����任��)��������ֵΪ��%f \n',value_max_change)

% [x111,fv111,exitflag1,outpt1,lama1] = linprog(-c1,A,b,[],[],zeros(length(c1),1));
% fv111=-fv111
% [x222,fv222,exitflag1,outpt1,lama1] = linprog(c2,A,b,[],[],zeros(length(c1),1));
% fv222
