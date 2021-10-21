clc
clear all
close all

%% 生成数据
tao=31;
N=40020;                  %设置点数N，Mackey_Glass默认时隙间隔为0.1;
[y,t]=Mackey_Glass(N,tao);%龙格-库塔(Runge-Kutta)方法求解Mackey Glass方程
n_train = 4000;
x_star = zeros(1,n_train);
for i = 1:n_train
    x_star(i) = y(1+i*10);
end

%% 相关的数据初始化

n = 4;%模糊规则前件个数为4
DATA_WIDE=1000;
DATA_SIZE=2;%注意此处只选1--2
Monte_carlo=20;
Tmax = 40;%此处不修改了40
RULE=[20,60];
RULE_MAX = max(RULE);
RULE_NUMBER = length(RULE);
CORE_NUMBER = 24;
POPULATION_NUM=4;
POPULATION_WIDE=40;

% n = 4;%模糊规则前件个数为4
% DATA_WIDE=500;
% DATA_SIZE=1;%注意此处只选1--2，太大费时间
% Monte_carlo=1;%此处设为1，不做蒙特卡罗实验了
% Tmax = 10;%此处不修改了
% RULE_WIDTH = 10;
% RULE_NUMBER = 1;
% CORE_NUMBER = 6;
% POPULATION_NUM=1;
% POPULATION_WIDE=10;
%% 遗传算法的相关参数
% NIND=20;          %种群大小
GA_MAXGEN=Tmax;        %最大遗传代数
GA_TIME=zeros(DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);
GA_parameter=zeros(3*(RULE_NUMBER*RULE_MAX)*n+2*(RULE_NUMBER*RULE_MAX),...
                   GA_MAXGEN,DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);
GA_RMSE=zeros(GA_MAXGEN,DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);
GA_RMSE_mean=zeros(GA_MAXGEN,DATA_SIZE,RULE_NUMBER,POPULATION_NUM);

%%

parpool(CORE_NUMBER);

%% 遗传算法（调用了谢菲尔德大学遗传算法工具箱的几个函数）

PRECI=20;         %个体长度
GGAP=0.95;        %代沟
px=0.8;           %交叉概率
pm=0.01;          %变异概率

for pop_num=1:POPULATION_NUM  %种群数的循环
    
    NIND=pop_num*POPULATION_WIDE;          %种群大小
    
for MC=1:Monte_carlo  %蒙特卡洛实验循环
    
    for rule_num=1:RULE_NUMBER   %规则数的循环
        
        M=RULE(rule_num);
        
        for i=1:DATA_SIZE        %数据规模的循环
           disp('***************************GA**********************************')
           fprintf(['种群数=',num2str(NIND),'蒙特卡洛次数=',num2str(MC),...
                    '规则数=',num2str(M),'数据规模=',num2str(i),'\n'])
        %% 选择数据
        X_train = [x_star(1:i*DATA_WIDE-4);x_star(2:i*DATA_WIDE-3);...
                   x_star(3:i*DATA_WIDE-2);x_star(4:i*DATA_WIDE-1);]';
        Y_train = x_star(5:i*DATA_WIDE)';
        X_test = [x_star(i*DATA_WIDE-3:(i+1)*DATA_WIDE-4);x_star(i*DATA_WIDE-2:(i+1)*DATA_WIDE-3);...
                  x_star(i*DATA_WIDE-1:(i+1)*DATA_WIDE-2);x_star(i*DATA_WIDE:(i+1)*DATA_WIDE-1)]';
        Y_test = x_star(i*DATA_WIDE+1:(i+1)*DATA_WIDE)';
        
        star_time=clock;  %开始计时
        
        trace=zeros(3*M*n+2*M+1,GA_MAXGEN);     %寻优结果的初始值
        FieldD=[PRECI*ones(1,3*M*n+2*M);rand(1,3*M*n+2*M);rand(1,3*M*n+2*M);...
            ones(1,3*M*n+2*M);zeros(1,3*M*n+2*M);ones(1,3*M*n+2*M);ones(1,3*M*n+2*M)];%区域描述器
        %M1= rand(M,n); M2= rand(M,n);sigma= rand(M,n);c1= rand(M,1); c2= rand(M,1);     
        Chrom = crtbp(NIND,PRECI*(3*M*n+2*M));  %创建任意离散随机种群
        %% 优化
        gen=0   %代计数器
        XY=bs2rv(Chrom,FieldD);  %初始种群的十进制转换
        ObjV=zeros(NIND,1);
        parfor nind_num=1:NIND   
               input=XY(nind_num,:);
               %[ M1,M2,sigma,c1,c2 ] = Mamdani_row_to_discrete( input,M,n );
               %sigma=sigma+0.00001;
               [R1,R2,R]=sfls_type2(X_train,input,M,n);
               ObjV(nind_num)=sum((Y_train-R').^2);%sqrt(sum((Y_train-R').^2)/length(R))
        end
        
        while gen<GA_MAXGEN
               if length(ObjV(:,1))==1   %ranking函数需要输入列向量，输入行向量报错，计算出NAN               
                  FitnV=ranking(ObjV');   %分配适应度值
               else
                  FitnV=ranking(ObjV);   %分配适应度值 
               end
              SelCh=select('sus',Chrom,FitnV,GGAP);  %选择
              SelCh=recombin('xovsp',SelCh,px);      %重组
              Selch=mut(SelCh,pm);                   %变异
              XY=bs2rv(SelCh,FieldD);                %子代个体的十进制转换
              ObjVSel=zeros(NIND,1);
              parfor nind_num=1:length(XY(:,1))
                     input=XY(nind_num,:);
                     %[ M1,M2,sigma,c1,c2 ] = Mamdani_row_to_discrete( input,M,n );
                     [R1,R2,R]=sfls_type2(X_train,input,M,n);
                     ObjVSel(nind_num)=sum((Y_train-R').^2);%sqrt(sum((Y_train-R').^2)/length(R))
              end   
              [Chrom,ObjV]=reins(Chrom,Selch,1,1,ObjV,ObjVSel(1:length(XY(:,1))));%重插入子代到父代，得到新种群
              XY=bs2rv(Chrom,FieldD);
              gen=gen+1
   
             %% 获取每一代的最优解及其序号，Y为最优解，I为个体的序号
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
        total_time=etime(end_time,star_time);  %结束计时
       %%
        GA_TIME(i,rule_num,MC,pop_num)=total_time;
        GA_P1=trace(1:3*M*n+2*M,:);
        GA_P2=trace(3*M*n+2*M+1,:);
        GA_parameter(1:length(GA_P1(:,1)),:,i,rule_num,MC,pop_num)=GA_P1;
        GA_RMSE(:,i,rule_num,MC,pop_num)= GA_P2;     
        end  %数据规模
    end %规则
    %delete(gcp('nocreate')) %每一轮循环后关闭MATLAB并行计算池，下一轮开启不同核心数的并行计算池
end %蒙特卡罗实验
 
   for MCMC=1:Monte_carlo
       GA_RMSE_mean(:,:,:,pop_num)=GA_RMSE_mean(:,:,:,pop_num)+GA_RMSE(:,:,:,MCMC,pop_num);
   end
   GA_RMSE_mean(:,:,:,pop_num)=GA_RMSE_mean(:,:,:,pop_num)./Monte_carlo;

end %种群数


delete(gcp('nocreate')) %关闭MATLAB并行计算池

% 保存生成的参数
save('GA_DATA','GA_RMSE_mean','GA_TIME','GA_parameter','GA_RMSE');














