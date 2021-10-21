function [xt_opt,Err,Error_train,time_train,t_awaystepFW] = ...
    Train_T1M_awaystepFW (xtrain,ytrain,x0,bound_low,bound_up,M,I,error_precision,epsilon,Tmax)
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
t_awaystepFW = 1;
Err = 20;
%T=1000;
Error_train = zeros(1,Tmax);
%rho = 0.5;
while (Err > error_precision) && t_awaystepFW <=Tmax    %ֹͣ����
    gama = 0.25/(t_awaystepFW+2);
    %gama = 0.1;
    if t_awaystepFW == 1
        xt = x0;
    end
    [ df, Err] = compute_df_f( xtrain,ytrain,xt,M,I,J );
   % Error_train = [Error_train Err];
    st = linprog(df,[],[],[],[],bound_low,bound_up);
    %st_save = [st_save st];
%% Away step FW
    d_FW=st-xt';
    vt = linprog(-df,[],[],[],[],bound_low,bound_up);
    d_A=xt'-vt;
    if -df*d_FW < epsilon
        break
    end
    if df*d_FW <= df*d_A
        dt=d_FW;%gama=1;
    else
        dt=d_A;%gama=0.5;
    end
%     gama_h=linspace(0.0001,gama,10);
%     for kh=1:length(gama_h)
%         [~,ERROR(kh)]=compute_df_f( xtrain(1:10,:),ytrain(1:10),xt+gama_h(kh)*dt',M,I,J );
%     end
%     locat=find(ERROR==min(ERROR));
%     xt=xt+gama_h(locat)*dt';
    
    xt=xt+gama*dt';
    [ ~, Err] = compute_df_f( xtrain,ytrain,xt,M,I,J );
    Err
    Error_train(t_awaystepFW)=Err;
    %xt_save = [xt_save xt'];
    %[ df1, Err1,fls_f1] = compute_df_f( xtrain,ytrain,xt,M,I,J );
    %Error_train = [Error_train Err1];
    t_awaystepFW = t_awaystepFW+1
    
    
%     figure(1)
%     plot([1:length(xtrain(:,1))],[ytrain fls_f1'],'LineWidth',1.0)
%     title(['��',num2str(t),'��ѵ��ʵ��ֵ��Ԥ��ֵͼ��'])
%     xlabel('t')
%     ylabel('y(t)')                    % ע��˴����������Լ�ͼƬ�������ڲ�ͬ��������Ҫ�޸�
%     legend('ʵ��ֵ','Ԥ��ֵ')  
   
end
%minerr_location = find(Error_train==min(Error_train)); %�ҵ�Error����С��λ��minerr_location
%xt_opt = xt_save(:,minerr_location(1))';
t_awaystepFW = t_awaystepFW -1;
xt_opt = xt;
time_train=toc;           %����ѵ��������ʱ��
end


