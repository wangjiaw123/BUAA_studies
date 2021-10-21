function [M1,M2,c1,c2,sigma,I2l,I2u,I1u,I1l]=train_sfls_type2(X,D,M1,M2,sigma,c1,c2,alpha)

[L,n]=size(X);
[N,n]=size(M1);

for i=1:L
U=[];
MU1=[];
UU=[];
LL=[];
for j=1:N
Uu=1;
Ll=1;
for m=1:n
P=[sigma(j,m),M1(j,m),M2(j,m)];
[uu,ll]=gausstype2(X(i,m),P);
Uu=Uu*uu;
Ll=Ll*ll;
end
UU=[UU,Uu];
LL=[LL,Ll];
end

%% Compute the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r_out,I2l,I2u,wr]= rightpoint(c2',LL,UU);          %此处用了EKM算法
[l_out,I1u,I1l,wl]= leftpoint(c1',LL,UU);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 此处用主单纯形法计算输出值
% n_psa = length(UU);
% A = zeros(2*n_psa,3*n_psa);
% A(1:n_psa,1:n_psa) = eye(n_psa);                     %此处用主单纯形法
% A(n_psa+1:2*n_psa,1:n_psa) = -eye(n_psa);
% A(1:2*n_psa,n_psa+1:3*n_psa) = eye(2*n_psa);
% b = [UU';-LL'];
% d1 = [ones(1,n_psa) zeros(1,2*n_psa)];
% d0 = 0;
% c0 = 0;
% c111 = [UU zeros(1,2*n_psa)]; %设置c1求最大值
% c222 = [LL zeros(1,2*n_psa)];%设置c2求最小值
% [ sheet1,sheet2,cow_location1,value_max_primal,cow_location_save1,T1,Xvalue1] = primal_simplex_algorithm( c111,c0,d1,d0,A,b );
% [ sheet11,sheet22,cow_location2,F2,cow_location_save2,T2,Xvalue2 ] = primal_simplex_algorithm( -c222,-c0,d1,d0,A,b );
% value_min_primal = -F2;
% r_out = value_max_primal;
% wr = Xvalue1(1:n_psa);
% I2l = [];I2u = [];
% for i = 1:n_psa
%     if abs(Xvalue1(i)-c111(i)) < abs(Xvalue1(i)-c222(i))
%         I2l = [I2l i];
%     elseif abs(Xvalue1(i)-c111(i)) >= abs(Xvalue1(i)-c222(i))
%         I2u = [I2u i];
%     end
% end
% 
% l_out = value_min_primal;
% wl = Xvalue2(1:n_psa);
% I1l = [];I1u = [];
% for i = 1:n_psa
%     if abs(Xvalue2(i)-c111(i)) <= abs(Xvalue2(i)-c222(i))
%         I1l = [I1l i];
%     elseif abs(Xvalue2(i)-c111(i)) > abs(Xvalue2(i)-c222(i))
%         I1u = [I1u i];
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 利用change_variable_method(变量变换法)计算输出值
% n_psa = length(UU);
% A = zeros(2*n_psa,3*n_psa);
% A(1:n_psa,1:n_psa) = eye(n_psa);                    
% A(n_psa+1:2*n_psa,1:n_psa) = -eye(n_psa);
% A(1:2*n_psa,n_psa+1:3*n_psa) = eye(2*n_psa);
% b = [UU';-LL'];
% d1 = [ones(1,n_psa) zeros(1,2*n_psa)];
% d0 = 0;
% c0 = 0;
% c111 = [UU zeros(1,2*n_psa)]; %设置c1求最大值
% c222 = [LL zeros(1,2*n_psa)];%设置c2求最小值
% AA = zeros(2*n_psa,n_psa);
% AA(1:n_psa,1:n_psa) = eye(n_psa);
% AA(n_psa+1:2*n_psa,1:n_psa) = -eye(n_psa);
% c11 = UU;
% d11 = ones(1,n_psa);
% % d11 = zeros(1,n_psa);
% % tic
% [ X_max_change,value_max_change ] = change_variable_method( c11,c0,d11,d0,AA,b );
% c22 = LL;
% [ X_min_change,value_min_change1 ] = change_variable_method( -c22,-c0,d11,d0,AA,b );
% % toc
% value_min_change = -value_min_change1;
% value_max_change;
% 
% r_out = value_max_change;
% wr = X_max_change(1:n_psa);
% I2l = [];I2u = [];
% for i = 1:n_psa
%     if abs(X_max_change(i)-c11(i)) < abs(X_max_change(i)-c22(i))
%         I2l = [I2l i];
%     elseif abs(X_max_change(i)-c11(i)) >= abs(X_max_change(i)-c22(i))
%         I2u = [I2u i];
%     end
% end
% 
% l_out = value_min_change;
% wl = X_min_change(1:n_psa);
% I1l = [];I1u = [];
% for i = 1:n_psa
%     if abs(X_min_change(i)-c11(i)) <= abs(X_min_change(i)-c22(i))
%         I1l = [I1l i];
%     elseif abs(X_min_change(i)-c11(i)) >= abs(X_min_change(i)-c22(i))
%         I1u = [I1u i];
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 利用iterative_method(迭代算法)计算输出值
% n_psa = length(UU);
% A = zeros(2*n_psa,3*n_psa);
% A(1:n_psa,1:n_psa) = eye(n_psa);                     
% A(n_psa+1:2*n_psa,1:n_psa) = -eye(n_psa);
% A(1:2*n_psa,n_psa+1:3*n_psa) = eye(2*n_psa);
% b = [UU';-LL'];
% d1 = [ones(1,n_psa) zeros(1,2*n_psa)];
% d0 = 0;
% c0 = 0;
% c111 = [UU zeros(1,2*n_psa)]; %设置c1求最大值
% c222 = [LL zeros(1,2*n_psa)];%设置c2求最小值
% tic
% [ X_max_iterative,value_max_iterative ] = iterative_method( c111,c0,d1,d0,A,b );
% [ X_min_iterative,value_min_iterative1 ] = iterative_method( -c222,-c0,d1,d0,A,b );
% toc
% value_min_iterative = -value_min_iterative1;
% r_out = value_max_iterative;
% wr = X_max_iterative(1:n_psa);
% I2l = [];I2u = [];
% for i = 1:n_psa
%     if abs(X_max_iterative(i)-c111(i)) < abs(X_max_iterative(i)-c222(i))
%         I2l = [I2l i];
%     elseif abs(X_max_iterative(i)-c111(i)) >= abs(X_max_iterative(i)-c222(i))
%         I2u = [I2u i];
%     end
% end
% 
% l_out = value_min_iterative;
% wl = X_min_iterative(1:n_psa);
% I1l = [];I1u = [];
% for i = 1:n_psa
%     if abs(X_min_iterative(i)-c111(i)) <= abs(X_min_iterative(i)-c222(i))
%         I1l = [I1l i];
%     elseif abs(X_min_iterative(i)-c111(i)) > abs(X_min_iterative(i)-c222(i))
%         I1u = [I1u i];
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f=(l_out+r_out)/2;
e=D(i)-f;

fa1=wl/sum(wl);
fa2=wr/sum(wr);

%%% Select MFs contributed to the right-most 
ME10=M1;
ME20=M2;
sigma0=sigma;

LE=length(I2u);

for t=1:LE
for k=1:n 
if X(i,k)<M1(I2u(t),k)
l=I2u(t);
M1(l,k)=M1(l,k)+alpha*e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(I2u(t))/sum(wr);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
elseif  X(i,k)>M2(I2u(t),k)
l=I2u(t);
M2(l,k)=M2(l,k)+alpha*e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(l)/sum(wr);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
end
end
end


LE=length(I2l);

for t=1:LE
for k=1:n
if X(i,k)<(M1(I2l(t),k)+M2(I2l(t),k))/2
l=I2l(t);
M2(l,k)=M2(l,k)+alpha*e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(l)/sum(wr);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
else
l=I2l(t);
M1(l,k)=M1(l,k)+alpha*e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c2(l)-r_out)*wr(l)/sum(wr);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c2(l)-r_out)*wr(l)/sum(wr);
end
end
end   

%%% Select MFs contributed to the left-most 

LE=length(I1l);

for t=1:LE
for k=1:n 
if X(i,k)<(M1(I1l(t),k)+M2(I1l(t),k))/2
l=I1l(t);
M2(l,k)=M2(l,k)+alpha*e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
else
l=I1l(t);
M1(l,k)=M1(l,k)+alpha*e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
end
end
end

LE=length(I1u);

for t=1:LE
for k=1:n
if  X(i,k)< M1(I1u(t),k)
l=I1u(t);
M1(l,k)=M1(l,k)+alpha*e*0.5*((X(i,k)-ME10(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME10(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
elseif  X(i,k)> M2(I1u(t),k)
l=I1u(t);
M2(l,k)=M2(l,k)+alpha*e*0.5*((X(i,k)-ME20(l,k))/(sigma(l,k)^2))...
*(c1(l)-l_out)*wl(l)/sum(wl);
sigma0(l,k)=sigma0(l,k)+alpha*e*0.5*(((X(i,k)-ME20(l,k))^2)/(sigma(l,k)^3))...
*(c1(l)-l_out)*wl(l)/sum(wl);
end
end
end  

sigma=sigma0;

for l=1:N
for k=1:n
if sigma(l,k) < 0
sigma(l,k)=abs(sigma(l,k));
end
end
end

fa1=wr'/sum(wr);
fa2=wl'/sum(wl);
c1=c1+alpha*e*fa1/2;
c2=c2+alpha*e*fa2/2;    


end

end