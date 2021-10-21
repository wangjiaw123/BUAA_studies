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
DATA_SIZE=2;%注意此处只选1--2，太大费时间6
Monte_carlo=20;%此处设为1，不做蒙特卡罗实验了
Tmax = 40;%此处不修改了
RULE=[20,60];
RULE_MAX = max(RULE);
RULE_NUMBER =  length(RULE);%4
CORE_NUMBER = 24;%24
POPULATION_NUM=4;%4
POPULATION_WIDE=40;%50

% % n = 4;%模糊规则前件个数为4
% DATA_WIDE=500;
% DATA_SIZE=1;%注意此处只选1--2，太大费时间
% Monte_carlo=1;%此处设为1，不做蒙特卡罗实验了
% Tmax = 10;%此处不修改了
% RULE_WIDTH = 10;
% RULE_NUMBER = 1;
% CORE_NUMBER = 6;
% POPULATION_NUM=1;
% POPULATION_WIDE=10;

%% 粒子群算法的相关参数
maxgen = Tmax; %迭代次数
%sizepop=20;  %种群规模
PS_TIME=zeros(DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);
PS_PARAMETER=zeros(Tmax,3*(RULE_NUMBER*RULE_MAX)*n+2*(RULE_NUMBER*RULE_MAX),...
                 DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);
PS_RMSE=zeros(Tmax,DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);             
PS_RMSE_MEAN=zeros(Tmax,DATA_SIZE,RULE_NUMBER,POPULATION_NUM);


parpool(CORE_NUMBER);

%% 粒子群算法(Particle swarm Algorithm)
%速度更新参数
pc1=1.49445;
pc2=1.49445;
%个体和速度最大最小值
popmax=1.5;popmin=0.001;
Vmax=0.1;Vmin=-0.1;

for pop_num=1:POPULATION_NUM %粒子数的循环
    
    sizepop=pop_num*POPULATION_WIDE; %种群数量规模
    
for MC=1:Monte_carlo  %蒙特卡洛实验循环
    
   
      
    for rule_num=1:RULE_NUMBER   %规则数的循环
        
        M=RULE(rule_num);
        
        for i=1:DATA_SIZE        %数据规模的循环
            disp('***************************PS**********************************')
            fprintf(['粒子数=',num2str(sizepop),',蒙特卡洛次数=',num2str(MC),...
                     ',规则数=',num2str(M),',数据规模=',num2str(i),'\n'])
        %% 选择数据
        X_train = [x_star(1:i*DATA_WIDE-4);x_star(2:i*DATA_WIDE-3);...
                   x_star(3:i*DATA_WIDE-2);x_star(4:i*DATA_WIDE-1);]';
        Y_train = x_star(5:i*DATA_WIDE)';
        X_test = [x_star(i*DATA_WIDE-3:(i+1)*DATA_WIDE-4);x_star(i*DATA_WIDE-2:(i+1)*DATA_WIDE-3);...
                  x_star(i*DATA_WIDE-1:(i+1)*DATA_WIDE-2);x_star(i*DATA_WIDE:(i+1)*DATA_WIDE-1)]';
        Y_test = x_star(i*DATA_WIDE+1:(i+1)*DATA_WIDE)';
        
        PS_star_time=clock;  %开始计时
       %% 种群初始化
        pop=zeros(sizepop,3*M*n+2*M);
        V=zeros(sizepop,3*M*n+2*M);
        fitness=zeros(1,sizepop);
        parfor ii=1:sizepop
             %随机产生一个种群 
             pop(ii,:)=rand(1,3*M*n+2*M);
             V(ii,:)=0.5*rand(1,3*M*n+2*M);
             input=pop(ii,:);
             [R1,R2,R]=sfls_type2(X_train,input,M,n);
             fitness(ii)=-sum((Y_train-R').^2);
        end
        %% 寻找初始极值（根据初始粒子适应度值寻找个体极值和群体极值）
        [bestfitness bestindex]=max(fitness);
        zbest=pop(bestindex(1),:);%群体极值位置
        gbest=pop;                %个体极值位置
        fitnessgbest=fitness;     %个体极值适应度值
        fitnesszbest=bestfitness; %群体极值适应度值
        %% 迭代寻优
        inside_RMSE=zeros(1,maxgen);
        inside_pare=zeros(maxgen,3*M*n+2*M);
        for ik=1:maxgen
            disp('//////////////PS//////////////')
            fprintf(['粒子数=',num2str(sizepop),',蒙特卡洛次数=',num2str(MC),...
                     ',规则数=',num2str(M),',数据规模=',num2str(i),'迭代次数=',num2str(ik),'\n'])
%%   此处代码这样写的话，parfor会将变量V分类错误                
%            for j=1:sizepop     %
%                j
%                %速度更新
%                V(j,:)=V(j,:)+pc1*rand(1)*(gbest(j,:)-pop(j,:))+pc2*rand(1)*(zbest-pop(j,:));
%                V(j,find(V(j,:)>Vmax))=Vmax;
%                V(j,find(V(j,:)<Vmin))=Vmin;
%                %粒子更新
%                pop(j,:)=pop(j,:)+0.5*V(j,:);
%                pop(j,find(pop(j,:)>popmax))=popmax;
%                pop(j,find(pop(j,:)<popmin))=popmin;               
%                %新粒子适应度值
%                disp('++++++++++++++++++++++++++')
%                input=pop(j,:);
%              %%  测试中，此处极易出错，算不出来，死循环               
%                [R1,R2,R]=sfls_type2(X_train,input+0.001*rand(size(input)),M,n);%
%                fitness(j)=-sum((Y_train-R').^2);
%                disp('+-+-+-+-+-+-+-+-+-+-+-+-+-')
%            end
%% 修改程序写法后，可以将此处并行
            V_help=zeros(size(V(1,:)));
            pop_help=zeros(size(pop(1,:)));
           parfor j=1:sizepop
               fprintf(['j=',num2str(j),'\n'])
               %% 分配
               V_help=V(j,:);
               pop_help=pop(j,:);
               %速度更新
               V_help=V_help+pc1*rand(1)*(gbest(j,:)-pop_help)+pc2*rand(1)*(zbest-pop_help);
               V_help(find(V_help>Vmax))=Vmax;
               V_help(find(V_help<Vmin))=Vmin;
               %粒子更新
               pop_help=pop_help+0.5*V_help;
               pop_help(find(pop_help>popmax))=popmax;
               pop_help(find(pop_help<popmin))=popmin;               
               %新粒子适应度值
               %disp('++++++++++++++++++++++++++')
               input=pop_help;
                %% 还原
               V(j,:)=V_help;
               pop(j,:)=pop_help;
             %%  测试中，此处极易出错，算不出来，死循环               
               [R1,R2,R]=sfls_type2(X_train,input+0.01*rand(size(input)),M,n);%
               fitness(j)=-sum((Y_train-R').^2);
               
               %disp('+-+-+-+-+-+-+-+-+-+-+-+-+-')
           end

               %% 个体极值和群体极值更新
               for j=1:sizepop
                   % 个体极值更新
                  if fitness(j)>fitnessgbest(j)
                      gbest(j,:)=pop(j,:);
                      fitnessgbest(j)=fitness(j);
                  end
                   % 群体极值更新
                  if fitness(j)>fitnesszbest
                      zbest=pop(j,:);
                      fitnesszbest=fitness(j);
                  end   
               end
               %disp('oooooooooooooooooooooooooooooo')
               % 记录每一代的最优值
               result(ik)=-fitnesszbest; 
              
               inside_pare(ik,:)=zbest;
               %% 记录预测时的RMSE
               [R1,R2,R]=sfls_type2_parfor(X_test,zbest,M,n);%+0.001*rand(size(input))             
               inside_RMSE(ik)=sqrt(sum(((Y_test-R').^2)/length(R)));              
        end
        
        PS_end_time=clock;
        PS_total_time=etime(PS_end_time,PS_star_time);  %结束计时
        
        PS_TIME(i,rule_num,MC,pop_num)=PS_total_time;
        PS_PARAMETER(:,1:length(inside_pare(1,:)),i,rule_num,MC,pop_num)=inside_pare;
        PS_RMSE(:,i,rule_num,MC,pop_num)=inside_RMSE;
        
        end %数据规模
   
    end %规则数

end %蒙特卡罗实验


for MCMC=1:Monte_carlo
   PS_RMSE_MEAN(:,:,:,pop_num) =PS_RMSE_MEAN(:,:,:,pop_num)+PS_RMSE(:,:,:,MCMC,pop_num); 
end
PS_RMSE_MEAN(:,:,:,pop_num)=PS_RMSE_MEAN(:,:,:,pop_num)/Monte_carlo;

end  %粒子数

delete(gcp('nocreate')) %关闭MATLAB并行计算池

% 保存粒子群算法产生的参数
save('PS_DATA','PS_RMSE_MEAN','PS_TIME','PS_RMSE','PS_PARAMETER')






