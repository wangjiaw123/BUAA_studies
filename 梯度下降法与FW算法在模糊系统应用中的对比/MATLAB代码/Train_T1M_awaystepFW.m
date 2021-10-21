function [xt_opt,Error_train,time_train,t_awaystepFW] = Train_T1M_awaystepFW (xtrain,ytrain,x0,bound_low,bound_up,M,I,error_precision,epsilon)
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
t_awaystepFW = 1;
Err = 200;
Error_train = [];
%rho = 0.5;
while Err > error_precision && t_awaystepFW <=5000    %停止条件
    gama = 2/(t_awaystepFW+2);
    if t_awaystepFW == 1
        xt = x0;
    end
    [ df, Err] = compute_df_f( xtrain,ytrain,xt,M,I,J );
   % Error_train = [Error_train Err];
    st = linprog(df,[],[],[],[],bound_low,bound_up);
    %st_save = [st_save st];
%%%%%%%%%%%%%%%%%%
%% 标准FW
%     gt = (df)*(st-xt');
%     xt = (1-gama)*xt+gama*st';
%     xt_save = [xt_save xt'];
%     %[ df1, Err1,fls_f1] = compute_df_f( xtrain,ytrain,xt,M,I,J );
%     %Error_train = [Error_train Err1];
%     if abs(norm(gt)) < epsilon
%         break;
%     end
%     t_FW = t_FW+1
%%%%%%%%%%%%%%%%%%
%% Away step FW
    d_FW=st-xt';
    vt = linprog(-df,[],[],[],[],bound_low,bound_up);
    d_A=xt'-vt;
    if -df*d_FW < epsilon
        %xt_opt=xt;
        break
    end
    if -df*d_FW >=-df*d_A
        dt=d_FW;
    else
        dt=d_A;
    end
    xt=xt+gama*dt';
    [ ~, Err] = compute_df_f( xtrain,ytrain,xt,M,I,J );
    %xt_save = [xt_save xt'];
    %[ df1, Err1,fls_f1] = compute_df_f( xtrain,ytrain,xt,M,I,J );
    %Error_train = [Error_train Err1];
    t_awaystepFW = t_awaystepFW+1
    
    
%     figure(1)
%     plot([1:length(xtrain(:,1))],[ytrain fls_f1'],'LineWidth',1.0)
%     title(['第',num2str(t),'次训练实际值与预测值图像'])
%     xlabel('t')
%     ylabel('y(t)')                    % 注意此处横纵坐标以及图片命名的在不同例子中需要修改
%     legend('实际值','预测值')  
   
end
%minerr_location = find(Error_train==min(Error_train)); %找到Error中最小的位置minerr_location
%xt_opt = xt_save(:,minerr_location(1))';
xt_opt = xt;
time_train=toc;           %计算训练参数的时间
end
