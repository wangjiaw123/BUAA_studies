function [xt_opt,Err,Error_train,time_train,t_FW] = ...
    Train_T1M_FW (xtrain,ytrain,x0,bound_low,bound_up,M,I,error_precision,epsilon,Tmax)
%%                   
%FW_FLS_train 使用FW算法训练一型单值模糊逻辑系统(FLS)的程序(规则为多前件一后
%              间，即令后件个数J=1)
%输入:  xtrain表示训练数据集(每个样本对应着行向量),ytrain表示训练数据集样本对
%       应的标签,M表示规则数目,I表示每条规则的前件个数，x0是初始化模糊系统的
%       参数，x0可看成行向量，其格式为[c(1),...,c(M),mu(1),...,mu(M*I),sigma(1),
%       ...,sigma(M*I)]。bound_low,bound_up是x0中参数对应的上下界,error_precision
%       表示跳出参数训练循化的最大误差，epsilon表示另一个停止条件。
%输出:  xt_opt表示最优的参数值,Error_train表示每次训练的误差，time_train表示训练时间
%%
J = 1;
tic        

%st_save = [];
%xt_save = [];
t_FW = 1;
Err = 20;
T=Tmax;%T=1000;
Error_train = zeros(1,Tmax);
%rho = 0.5;
while Err > error_precision  && t_FW<=T   %停止条件
    gama = 2/(t_FW+2);
    if t_FW == 1
        xt = x0;
    end
    [ df, Err] = compute_df_f( xtrain,ytrain,xt,M,I,J );
   % Error_train = [Error_train Err];
    st = linprog(df,[],[],[],[],bound_low,bound_up);

%% 标准FW
    gt = (df)*(st-xt');
    xt = (1-gama)*xt+gama*st';
    %xt_save = [xt_save xt'];
    [ ~, Err,~] = compute_df_f( xtrain,ytrain,xt,M,I,J );
    Error_train(t_FW)=Err;
    %Error_train = [Error_train Err1];
    if abs(norm(gt)) < epsilon
        break;
    end
    t_FW = t_FW+1   
end
%minerr_location = find(Error_train==min(Error_train)); %找到Error中最小的位置minerr_location
%xt_opt = xt_save(:,minerr_location(1))';
t_FW = t_FW-1;
xt_opt = xt;
time_train=toc;           %计算训练参数的时间
end

