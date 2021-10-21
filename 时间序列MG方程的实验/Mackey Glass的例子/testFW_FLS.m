close all
clear
clc
%% 生成数据
tao=31;
Dt=tao;
N=20020;%设置点数N，时间间隔为0.1;
[y,t]=Mackey_Glass(N,tao);%龙格-库塔(Runge-Kutta)方法求解Mackey Glass方程
% figure(1)
% subplot(1,2,1)% 序列直出 
% plot(t,y,'LineWidth',1.0);
% xlabel('t')
% ylabel('s(t)')
% title('序列直出 ');
% subplot(1,2,2)% 相图时差 
% plot(y((Dt*10+1+10000):end),y(10001:end-10*Dt),'LineWidth',0.5);
% xlabel('s(t)')
% ylabel('s(t-\tau)')
% title('相图时差');

n_tr = 1000;
x_star = zeros(1,n_tr);
for i = 1:n_tr
    x_star(i) = y(10001+i*10);
end
%% 训练数据以及测试数据对
X_train = [x_star(1:500);x_star(2:501);x_star(3:502);x_star(4:503);]';
Y_train = x_star(5:504)';
X_test = [x_star(505:996);x_star(506:997);x_star(507:998);x_star(508:999);]';
Y_test = x_star(509:1000)';
%% 
n = length(X_train(1,:));
N = 20;%设置规则数目20
M1 = rand(N,n);
M2 = M1 + 2;0.3
sigma = rand(N,n);
c1 = rand(N,1);
c2 = c1+0.55;0.6/0.55
cc1 = c1;
cc2 = c2;

rho=0.4;%%%%%%%
alpha=0.5;%%%%% 两个是Armijo非线性线搜索的参数
epsilon=10e-3;
st_save=[];
Error=[];
x0=rand(1,3*N*n+2*N);
% M1re1=rand(3*n*N);
% M2re2=M1re1+2; 
% sigma1re1=10*rand(1,3*n*N);
% c1re1=0.5*rand(1,N);
% c2re2=clre1+2;
%bound_low=[0.1*ones(1,n*N),0.5*ones(1,n*N),0.5*ones(1,n*N),-1*rand(1,N),1*ones(1,N)];
bound_low=[0.1*rand(1,n*N),0.1*rand(1,n*N),0.1*ones(1,n*N),0.1*rand(1,N),0.5*ones(1,N)];
bound_up=[1*ones(1,n*N),1*ones(1,n*N),4*ones(1,n*N),0.4*ones(1,N),1*ones(1,N)];
% bound_low=[0.1*ones(1,n*N),0.5*ones(1,n*N),0.5*ones(1,n*N),0.1*ones(1,N),0.1*ones(1,N)];
% bound_up=[3*rand(1,n*N),3*rand(1,n*N),3*rand(1,n*N),0.5*ones(1,N)+1,0.5*ones(1,N)+1];
% bound_low=0.1*ones(1,3*N*n+2*N);%rand(3*n*N)
% bound_up=bound_low+0.5;
Tmax=12;
t=1;
while t<Tmax %限制迭代次数为**
    if t==1
        xt=x0;
    end
%% 计算总的误差Err
%    [ df, Err,fls_f] = compute_df_f( xtrain,ytrain,xt,M,I,J );
    [Err]=FW_FLS_Err(X_train,Y_train,M1,M2,sigma,c1,c2)
    Error=[Error Err];
 %% 计算误差对各个变量的偏导数
    M111=zeros(size(M1));
    M222=zeros(size(M2));
    sigma111=zeros(size(sigma));
    c111=zeros(size(c1));
    c222=zeros(size(c2));
    df=[];
    for i=1:length(Y_train)
        [M11,M22,c11,c22,sigma11]=fun_derivation(X_train(i,:),Y_train(i),M1,M2,sigma,c1,c2);
        M111=M111+M11;
        M222=M222+M22;
        sigma111=sigma111+sigma11;
        c111=c111+c11;
        c222=c222+c22;

    end
    M1re=[];
    M2re=[];
    sigmalre=[];
    c1re=[];
    c2re=[];
    M1re=reshape(M111',1,[]);
    M2re=reshape(M222',1,[]); 
    sigma1re=reshape(sigma11',1,[]);
    c1re=reshape(c111,1,[]);
    c2re=reshape(c222,1,[]);
    df=[M1re,M2re,sigma1re,c1re,c2re];
%%    
    st=linprog(df,[],[],[],[],bound_low,bound_up);
    st_save(:,t)=st;
    dt = st-xt';
    gt = (df)*dt;
    if abs(gt) < epsilon
        break;
    end
    m=0;
    while m< 20 %非精确线搜索的Armijo方法
        [M1,M2,sigma,c1,c2] = transform(xt,n,N);
        [Err]=FW_FLS_Err(X_train,Y_train,M1,M2,sigma,c1,c2);
%        [df2,Err]=compute_df_f( xtrain,ytrain,xt,M,I,J );
        fvalue0=Err; 
        xt=xt+rho^m*dt';
        [M17,M27,sigma17,c17,c27] = transform(xt,n,N);
        [Err]=FW_FLS_Err(X_train,Y_train,M17,M27,sigma17,c17,c27);
%        [df1,Err]=compute_df_f( X_train,Y_train,xt,M,I,J );
        fvalue1=Err;
        mk=0;
        if fvalue1<fvalue0+alpha*rho^m*df*dt%(-df')
             mk=m; 
             break;
        end
        m=m+1;
    end
    xt=xt+rho^mk*dt'; 
    [M17,M27,sigma17,c17,c27] = transform(xt,n,N);
    [R1,R2,R]=sfls_type2(X_train,M1,M2,sigma,c1,c2);
    figure(2)
    plot(1:length(X_train),[Y_train';R]);
    title(['第',num2str(t),'次训练实际值与预测值图像','(\tau =',num2str(tao),')'])
    xlabel('t')
    ylabel('s(t)')
    legend('实际值','预测值')
    t=t+1    
%     c1=c1';c2=c2';
end   

figure
plot(Error)
legend('损失(目标)函数')
minerr_location=find(Error==min(Error));%找到Error中最小的位置minerr_location

[M1,M2,sigma,c1,c2] = transform(st_save(:,minerr_location),n,N);
[R1,R2,R]=sfls_type2(X_test,M1,M2,sigma,c1',c2');
figure 
plot(1500+1:1500+length(X_test),R)
hold on
plot(1500+1:1500+length(X_test),Y_test)
xlabel('t')
ylabel('s(t)')
legend('预测值','实际值')             
% plot(fls_f)
% hold on
% plot(y)
% legend('预测','实际图像')
% title('用F-W算法优化模糊系统参数后的预测输出与实际图像')
% figure
% plot(fls_f)
% legend('预测')
%sprintf('%f',Error(minerr_location))
strvcat('平方误差',num2str(Error(minerr_location)))
disp(sprintf('平方误差为%d',Error(minerr_location)))




























