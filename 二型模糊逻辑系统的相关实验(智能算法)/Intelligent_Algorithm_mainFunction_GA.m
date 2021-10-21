clc
clear all
close all

%% ��������
tao=31;
N=40020;                  %���õ���N��Mackey_GlassĬ��ʱ϶���Ϊ0.1;
[y,t]=Mackey_Glass(N,tao);%����-����(Runge-Kutta)�������Mackey Glass����
n_train = 4000;
x_star = zeros(1,n_train);
for i = 1:n_train
    x_star(i) = y(1+i*10);
end

%% ��ص����ݳ�ʼ��

n = 4;%ģ������ǰ������Ϊ4
DATA_WIDE=1000;
DATA_SIZE=2;%ע��˴�ֻѡ1--2
Monte_carlo=20;
Tmax = 40;%�˴����޸���40
RULE=[20,60];
RULE_MAX = max(RULE);
RULE_NUMBER = length(RULE);
CORE_NUMBER = 24;
POPULATION_NUM=4;
POPULATION_WIDE=40;

% n = 4;%ģ������ǰ������Ϊ4
% DATA_WIDE=500;
% DATA_SIZE=1;%ע��˴�ֻѡ1--2��̫���ʱ��
% Monte_carlo=1;%�˴���Ϊ1���������ؿ���ʵ����
% Tmax = 10;%�˴����޸���
% RULE_WIDTH = 10;
% RULE_NUMBER = 1;
% CORE_NUMBER = 6;
% POPULATION_NUM=1;
% POPULATION_WIDE=10;
%% �Ŵ��㷨����ز���
% NIND=20;          %��Ⱥ��С
GA_MAXGEN=Tmax;        %����Ŵ�����
GA_TIME=zeros(DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);
GA_parameter=zeros(3*(RULE_NUMBER*RULE_MAX)*n+2*(RULE_NUMBER*RULE_MAX),...
                   GA_MAXGEN,DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);
GA_RMSE=zeros(GA_MAXGEN,DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);
GA_RMSE_mean=zeros(GA_MAXGEN,DATA_SIZE,RULE_NUMBER,POPULATION_NUM);

%%

parpool(CORE_NUMBER);

%% �Ŵ��㷨��������л�ƶ��´�ѧ�Ŵ��㷨������ļ���������

PRECI=20;         %���峤��
GGAP=0.95;        %����
px=0.8;           %�������
pm=0.01;          %�������

for pop_num=1:POPULATION_NUM  %��Ⱥ����ѭ��
    
    NIND=pop_num*POPULATION_WIDE;          %��Ⱥ��С
    
for MC=1:Monte_carlo  %���ؿ���ʵ��ѭ��
    
    for rule_num=1:RULE_NUMBER   %��������ѭ��
        
        M=RULE(rule_num);
        
        for i=1:DATA_SIZE        %���ݹ�ģ��ѭ��
           disp('***************************GA**********************************')
           fprintf(['��Ⱥ��=',num2str(NIND),'���ؿ������=',num2str(MC),...
                    '������=',num2str(M),'���ݹ�ģ=',num2str(i),'\n'])
        %% ѡ������
        X_train = [x_star(1:i*DATA_WIDE-4);x_star(2:i*DATA_WIDE-3);...
                   x_star(3:i*DATA_WIDE-2);x_star(4:i*DATA_WIDE-1);]';
        Y_train = x_star(5:i*DATA_WIDE)';
        X_test = [x_star(i*DATA_WIDE-3:(i+1)*DATA_WIDE-4);x_star(i*DATA_WIDE-2:(i+1)*DATA_WIDE-3);...
                  x_star(i*DATA_WIDE-1:(i+1)*DATA_WIDE-2);x_star(i*DATA_WIDE:(i+1)*DATA_WIDE-1)]';
        Y_test = x_star(i*DATA_WIDE+1:(i+1)*DATA_WIDE)';
        
        star_time=clock;  %��ʼ��ʱ
        
        trace=zeros(3*M*n+2*M+1,GA_MAXGEN);     %Ѱ�Ž���ĳ�ʼֵ
        FieldD=[PRECI*ones(1,3*M*n+2*M);rand(1,3*M*n+2*M);rand(1,3*M*n+2*M);...
            ones(1,3*M*n+2*M);zeros(1,3*M*n+2*M);ones(1,3*M*n+2*M);ones(1,3*M*n+2*M)];%����������
        %M1= rand(M,n); M2= rand(M,n);sigma= rand(M,n);c1= rand(M,1); c2= rand(M,1);     
        Chrom = crtbp(NIND,PRECI*(3*M*n+2*M));  %����������ɢ�����Ⱥ
        %% �Ż�
        gen=0   %��������
        XY=bs2rv(Chrom,FieldD);  %��ʼ��Ⱥ��ʮ����ת��
        ObjV=zeros(NIND,1);
        parfor nind_num=1:NIND   
               input=XY(nind_num,:);
               %[ M1,M2,sigma,c1,c2 ] = Mamdani_row_to_discrete( input,M,n );
               %sigma=sigma+0.00001;
               [R1,R2,R]=sfls_type2(X_train,input,M,n);
               ObjV(nind_num)=sum((Y_train-R').^2);%sqrt(sum((Y_train-R').^2)/length(R))
        end
        
        while gen<GA_MAXGEN
               if length(ObjV(:,1))==1   %ranking������Ҫ�������������������������������NAN               
                  FitnV=ranking(ObjV');   %������Ӧ��ֵ
               else
                  FitnV=ranking(ObjV);   %������Ӧ��ֵ 
               end
              SelCh=select('sus',Chrom,FitnV,GGAP);  %ѡ��
              SelCh=recombin('xovsp',SelCh,px);      %����
              Selch=mut(SelCh,pm);                   %����
              XY=bs2rv(SelCh,FieldD);                %�Ӵ������ʮ����ת��
              ObjVSel=zeros(NIND,1);
              parfor nind_num=1:length(XY(:,1))
                     input=XY(nind_num,:);
                     %[ M1,M2,sigma,c1,c2 ] = Mamdani_row_to_discrete( input,M,n );
                     [R1,R2,R]=sfls_type2(X_train,input,M,n);
                     ObjVSel(nind_num)=sum((Y_train-R').^2);%sqrt(sum((Y_train-R').^2)/length(R))
              end   
              [Chrom,ObjV]=reins(Chrom,Selch,1,1,ObjV,ObjVSel(1:length(XY(:,1))));%�ز����Ӵ����������õ�����Ⱥ
              XY=bs2rv(Chrom,FieldD);
              gen=gen+1
   
             %% ��ȡÿһ�������Ž⼰����ţ�YΪ���Ž⣬IΪ��������
              %[Y,I]=max(ObjV);
              [Y,I]=min(ObjV);
              
              predit_ObjVSel=zeros(NIND,1);

              input=XY(I(1),:);
              [ M1,M2,sigma,c1,c2 ] = Mamdani_row_to_discrete( input,M,n );
              sigma=sigma+0.00001;
              [R1,R2,R]=sfls_type2_parfor(X_test,input,M,n);
              RMSE= sqrt(sum((Y_test-R').^2)/length(R))
              trace(1:3*M*n+2*M,gen)=XY(I(1),:);
              trace(3*M*n+2*M+1,gen)=RMSE;  

        end
        end_time=clock;
        total_time=etime(end_time,star_time);  %������ʱ
       %%
        GA_TIME(i,rule_num,MC,pop_num)=total_time;
        GA_P1=trace(1:3*M*n+2*M,:);
        GA_P2=trace(3*M*n+2*M+1,:);
        GA_parameter(1:length(GA_P1(:,1)),:,i,rule_num,MC,pop_num)=GA_P1;
        GA_RMSE(:,i,rule_num,MC,pop_num)= GA_P2;     
        end  %���ݹ�ģ
    end %����
    %delete(gcp('nocreate')) %ÿһ��ѭ����ر�MATLAB���м���أ���һ�ֿ�����ͬ�������Ĳ��м����
end %���ؿ���ʵ��
 
   for MCMC=1:Monte_carlo
       GA_RMSE_mean(:,:,:,pop_num)=GA_RMSE_mean(:,:,:,pop_num)+GA_RMSE(:,:,:,MCMC,pop_num);
   end
   GA_RMSE_mean(:,:,:,pop_num)=GA_RMSE_mean(:,:,:,pop_num)./Monte_carlo;

end %��Ⱥ��


delete(gcp('nocreate')) %�ر�MATLAB���м����

% �������ɵĲ���
save('GA_DATA','GA_RMSE_mean','GA_TIME','GA_parameter','GA_RMSE');














