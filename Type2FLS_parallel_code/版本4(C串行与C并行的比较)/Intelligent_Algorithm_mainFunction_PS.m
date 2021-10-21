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
DATA_SIZE=2;%ע��˴�ֻѡ1--2��̫���ʱ��6
Monte_carlo=20;%�˴���Ϊ1���������ؿ���ʵ����
Tmax = 40;%�˴����޸���
RULE=[20,60];
RULE_MAX = max(RULE);
RULE_NUMBER =  length(RULE);%4
CORE_NUMBER = 24;%24
POPULATION_NUM=4;%4
POPULATION_WIDE=40;%50

% % n = 4;%ģ������ǰ������Ϊ4
% DATA_WIDE=500;
% DATA_SIZE=1;%ע��˴�ֻѡ1--2��̫���ʱ��
% Monte_carlo=1;%�˴���Ϊ1���������ؿ���ʵ����
% Tmax = 10;%�˴����޸���
% RULE_WIDTH = 10;
% RULE_NUMBER = 1;
% CORE_NUMBER = 6;
% POPULATION_NUM=1;
% POPULATION_WIDE=10;

%% ����Ⱥ�㷨����ز���
maxgen = Tmax; %��������
%sizepop=20;  %��Ⱥ��ģ
PS_TIME=zeros(DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);
PS_PARAMETER=zeros(Tmax,3*(RULE_NUMBER*RULE_MAX)*n+2*(RULE_NUMBER*RULE_MAX),...
                 DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);
PS_RMSE=zeros(Tmax,DATA_SIZE,RULE_NUMBER,Monte_carlo,POPULATION_NUM);             
PS_RMSE_MEAN=zeros(Tmax,DATA_SIZE,RULE_NUMBER,POPULATION_NUM);


parpool(CORE_NUMBER);

%% ����Ⱥ�㷨(Particle swarm Algorithm)
%�ٶȸ��²���
pc1=1.49445;
pc2=1.49445;
%������ٶ������Сֵ
popmax=1.5;popmin=0.001;
Vmax=0.1;Vmin=-0.1;

for pop_num=1:POPULATION_NUM %��������ѭ��
    
    sizepop=pop_num*POPULATION_WIDE; %��Ⱥ������ģ
    
for MC=1:Monte_carlo  %���ؿ���ʵ��ѭ��
    
   
      
    for rule_num=1:RULE_NUMBER   %��������ѭ��
        
        M=RULE(rule_num);
        
        for i=1:DATA_SIZE        %���ݹ�ģ��ѭ��
            disp('***************************PS**********************************')
            fprintf(['������=',num2str(sizepop),',���ؿ������=',num2str(MC),...
                     ',������=',num2str(M),',���ݹ�ģ=',num2str(i),'\n'])
        %% ѡ������
        X_train = [x_star(1:i*DATA_WIDE-4);x_star(2:i*DATA_WIDE-3);...
                   x_star(3:i*DATA_WIDE-2);x_star(4:i*DATA_WIDE-1);]';
        Y_train = x_star(5:i*DATA_WIDE)';
        X_test = [x_star(i*DATA_WIDE-3:(i+1)*DATA_WIDE-4);x_star(i*DATA_WIDE-2:(i+1)*DATA_WIDE-3);...
                  x_star(i*DATA_WIDE-1:(i+1)*DATA_WIDE-2);x_star(i*DATA_WIDE:(i+1)*DATA_WIDE-1)]';
        Y_test = x_star(i*DATA_WIDE+1:(i+1)*DATA_WIDE)';
        
        PS_star_time=clock;  %��ʼ��ʱ
       %% ��Ⱥ��ʼ��
        pop=zeros(sizepop,3*M*n+2*M);
        V=zeros(sizepop,3*M*n+2*M);
        fitness=zeros(1,sizepop);
        parfor ii=1:sizepop
             %�������һ����Ⱥ 
             pop(ii,:)=rand(1,3*M*n+2*M);
             V(ii,:)=0.5*rand(1,3*M*n+2*M);
             input=pop(ii,:);
             [R1,R2,R]=sfls_type2(X_train,input,M,n);
             fitness(ii)=-sum((Y_train-R').^2);
        end
        %% Ѱ�ҳ�ʼ��ֵ�����ݳ�ʼ������Ӧ��ֵѰ�Ҹ��弫ֵ��Ⱥ�弫ֵ��
        [bestfitness bestindex]=max(fitness);
        zbest=pop(bestindex(1),:);%Ⱥ�弫ֵλ��
        gbest=pop;                %���弫ֵλ��
        fitnessgbest=fitness;     %���弫ֵ��Ӧ��ֵ
        fitnesszbest=bestfitness; %Ⱥ�弫ֵ��Ӧ��ֵ
        %% ����Ѱ��
        inside_RMSE=zeros(1,maxgen);
        inside_pare=zeros(maxgen,3*M*n+2*M);
        for ik=1:maxgen
            disp('//////////////PS//////////////')
            fprintf(['������=',num2str(sizepop),',���ؿ������=',num2str(MC),...
                     ',������=',num2str(M),',���ݹ�ģ=',num2str(i),'��������=',num2str(ik),'\n'])
%%   �˴���������д�Ļ���parfor�Ὣ����V�������                
%            for j=1:sizepop     %
%                j
%                %�ٶȸ���
%                V(j,:)=V(j,:)+pc1*rand(1)*(gbest(j,:)-pop(j,:))+pc2*rand(1)*(zbest-pop(j,:));
%                V(j,find(V(j,:)>Vmax))=Vmax;
%                V(j,find(V(j,:)<Vmin))=Vmin;
%                %���Ӹ���
%                pop(j,:)=pop(j,:)+0.5*V(j,:);
%                pop(j,find(pop(j,:)>popmax))=popmax;
%                pop(j,find(pop(j,:)<popmin))=popmin;               
%                %��������Ӧ��ֵ
%                disp('++++++++++++++++++++++++++')
%                input=pop(j,:);
%              %%  �����У��˴����׳����㲻��������ѭ��               
%                [R1,R2,R]=sfls_type2(X_train,input+0.001*rand(size(input)),M,n);%
%                fitness(j)=-sum((Y_train-R').^2);
%                disp('+-+-+-+-+-+-+-+-+-+-+-+-+-')
%            end
%% �޸ĳ���д���󣬿��Խ��˴�����
            V_help=zeros(size(V(1,:)));
            pop_help=zeros(size(pop(1,:)));
           parfor j=1:sizepop
               fprintf(['j=',num2str(j),'\n'])
               %% ����
               V_help=V(j,:);
               pop_help=pop(j,:);
               %�ٶȸ���
               V_help=V_help+pc1*rand(1)*(gbest(j,:)-pop_help)+pc2*rand(1)*(zbest-pop_help);
               V_help(find(V_help>Vmax))=Vmax;
               V_help(find(V_help<Vmin))=Vmin;
               %���Ӹ���
               pop_help=pop_help+0.5*V_help;
               pop_help(find(pop_help>popmax))=popmax;
               pop_help(find(pop_help<popmin))=popmin;               
               %��������Ӧ��ֵ
               %disp('++++++++++++++++++++++++++')
               input=pop_help;
                %% ��ԭ
               V(j,:)=V_help;
               pop(j,:)=pop_help;
             %%  �����У��˴����׳����㲻��������ѭ��               
               [R1,R2,R]=sfls_type2(X_train,input+0.01*rand(size(input)),M,n);%
               fitness(j)=-sum((Y_train-R').^2);
               
               %disp('+-+-+-+-+-+-+-+-+-+-+-+-+-')
           end

               %% ���弫ֵ��Ⱥ�弫ֵ����
               for j=1:sizepop
                   % ���弫ֵ����
                  if fitness(j)>fitnessgbest(j)
                      gbest(j,:)=pop(j,:);
                      fitnessgbest(j)=fitness(j);
                  end
                   % Ⱥ�弫ֵ����
                  if fitness(j)>fitnesszbest
                      zbest=pop(j,:);
                      fitnesszbest=fitness(j);
                  end   
               end
               %disp('oooooooooooooooooooooooooooooo')
               % ��¼ÿһ��������ֵ
               result(ik)=-fitnesszbest; 
              
               inside_pare(ik,:)=zbest;
               %% ��¼Ԥ��ʱ��RMSE
               [R1,R2,R]=sfls_type2_parfor(X_test,zbest,M,n);%+0.001*rand(size(input))             
               inside_RMSE(ik)=sqrt(sum(((Y_test-R').^2)/length(R)));              
        end
        
        PS_end_time=clock;
        PS_total_time=etime(PS_end_time,PS_star_time);  %������ʱ
        
        PS_TIME(i,rule_num,MC,pop_num)=PS_total_time;
        PS_PARAMETER(:,1:length(inside_pare(1,:)),i,rule_num,MC,pop_num)=inside_pare;
        PS_RMSE(:,i,rule_num,MC,pop_num)=inside_RMSE;
        
        end %���ݹ�ģ
   
    end %������

end %���ؿ���ʵ��


for MCMC=1:Monte_carlo
   PS_RMSE_MEAN(:,:,:,pop_num) =PS_RMSE_MEAN(:,:,:,pop_num)+PS_RMSE(:,:,:,MCMC,pop_num); 
end
PS_RMSE_MEAN(:,:,:,pop_num)=PS_RMSE_MEAN(:,:,:,pop_num)/Monte_carlo;

end  %������

delete(gcp('nocreate')) %�ر�MATLAB���м����

% ��������Ⱥ�㷨�����Ĳ���
save('PS_DATA','PS_RMSE_MEAN','PS_TIME','PS_RMSE','PS_PARAMETER')






