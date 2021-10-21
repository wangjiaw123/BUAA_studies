function [xt_opt,Err,Error_train,time_train,t_FW] = ...
    Train_T1M_FW (xtrain,ytrain,x0,bound_low,bound_up,M,I,error_precision,epsilon,Tmax)
%%                   
%FW_FLS_train ʹ��FW�㷨ѵ��һ�͵�ֵģ���߼�ϵͳ(FLS)�ĳ���(����Ϊ��ǰ��һ��
%              �䣬����������J=1)
%����:  xtrain��ʾѵ�����ݼ�(ÿ��������Ӧ��������),ytrain��ʾѵ�����ݼ�������
%       Ӧ�ı�ǩ,M��ʾ������Ŀ,I��ʾÿ�������ǰ��������x0�ǳ�ʼ��ģ��ϵͳ��
%       ������x0�ɿ��������������ʽΪ[c(1),...,c(M),mu(1),...,mu(M*I),sigma(1),
%       ...,sigma(M*I)]��bound_low,bound_up��x0�в�����Ӧ�����½�,error_precision
%       ��ʾ��������ѵ��ѭ���������epsilon��ʾ��һ��ֹͣ������
%���:  xt_opt��ʾ���ŵĲ���ֵ,Error_train��ʾÿ��ѵ������time_train��ʾѵ��ʱ��
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
while Err > error_precision  && t_FW<=T   %ֹͣ����
    gama = 2/(t_FW+2);
    if t_FW == 1
        xt = x0;
    end
    [ df, Err] = compute_df_f( xtrain,ytrain,xt,M,I,J );
   % Error_train = [Error_train Err];
    st = linprog(df,[],[],[],[],bound_low,bound_up);

%% ��׼FW
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
%minerr_location = find(Error_train==min(Error_train)); %�ҵ�Error����С��λ��minerr_location
%xt_opt = xt_save(:,minerr_location(1))';
t_FW = t_FW-1;
xt_opt = xt;
time_train=toc;           %����ѵ��������ʱ��
end

