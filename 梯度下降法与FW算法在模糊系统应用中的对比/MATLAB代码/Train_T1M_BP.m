function [meanF,stdF,yl,V,S]=Train_T1M_BP(input,output,meanF,stdF,yl,alpha,V,S)

% The function T1TSK is to get parameters of a TSK-M FLS from given inputs,outputs
% input: training inputs
% output: training outputs
% meanF: mean of antecedent mfs F_k^l
% stdF: standard deviation of antecedent mfs F_k^l
% yl: centers of consequent fuzzy sets
% alpha: learning rate

[T,~]=size(input);
[M,n]=size(meanF);

%% Inference by using Mamdani FLS and tuning parameters by BP 

fl=zeros(M,1); % Firing levels f^l
output1=zeros(T,1);
%% 梯度下降法
    for t=1:T
        % Inference by using TSK FLS
        for l=1:M
            u=1;
            for i=1:n
                u1=exp(-(input(t,i)-meanF(l,i))^2/(2*stdF(l,i)^2));
                u=u*u1;
            end
            fl(l)=u;
        end
        phi=fl/sum(fl);
        output1(t)=yl'*phi; 
    
        % Tuning parameters by BP
    
        meanF1=zeros(M,n);
        stdF1=zeros(M,n);
        yl1=yl-alpha*(output1(t)-output(t))*phi;
        
        for i=1:n
            meanF1(:,i)=meanF(:,i)-alpha*(output1(t)-output(t))*(yl-output1(t)).*phi.*((input(t,i)-meanF(:,i))./stdF(:,i).^2);
            stdF1(:,i)=stdF(:,i)-alpha*(output1(t)-output(t)).*(yl-output1(t)).*phi.*(input(t,i)-meanF(:,i)).^2./stdF(:,i).^3;
        end
        yl=yl1;meanF=meanF1;stdF=stdF1;
    end
    
%% 随机梯度下降法(速度太慢，收敛不了)    
%     t=randi([1,T],1,1);
%         % Inference by using TSK FLS
%         for l=1:M
%             u=1;
%             for i=1:n
%                 u1=exp(-(input(t,i)-meanF(l,i))^2/(2*stdF(l,i)^2));
%                 u=u*u1;
%             end
%             fl(l)=u;
%         end
%         phi=fl/sum(fl);
%         output1(t)=yl'*phi; 
%     
%         % Tuning parameters by BP
%     
%         meanF1=zeros(M,n);
%         stdF1=zeros(M,n);
%         yl=yl-alpha*(output1(t)-output(t))*phi;
%         meanF=meanF-alpha*(output1(t)-output(t))*(yl-output1(t)).*phi.*((input(t,i)-meanF(:,i))./stdF(:,i).^2);
%         stdF=stdF-alpha*(output1(t)-output(t)).*(yl-output1(t)).*phi.*(input(t,i)-meanF(:,i)).^2./stdF(:,i).^3;
%         %yl=yl1;meanF=meanF1;stdF=stdF1;
%%  Adam算法(有错误)
% xt=[yl',reshape(meanF',[1,M*n]),reshape(stdF',[1,M*n])];
% beta1=0.9;
% beta2=0.999
%     for t=1:T
%         if t==1
%             V=zeros(size(1,(2*n+1)*M));
%             S=zeros(size(1,(2*n+1)*M));
%         end
%         % Inference by using TSK FLS
%         for l=1:M
%             u=1;
%             for i=1:n
%                 u1=exp(-(input(t,i)-meanF(l,i))^2/(2*stdF(l,i)^2));
%                 u=u*u1;
%             end
%             fl(l)=u;
%         end
%         phi=fl/sum(fl);
%         output1(t)=yl'*phi; 
%     
%         % Tuning parameters by BP
%     
%         meanF1=zeros(M,n);
%         stdF1=zeros(M,n);
%         yll=zeros(size(yl));
%         meanF1=zeros(size(meanF));
%         stdF1=zeros(size(stdF));
%         
%         yll=yll+(output1(t)-output(t))*phi;
%         for i=1:n
%             meanF1(:,i)=meanF1(:,i)+(output1(t)-output(t))*(yl-output1(t)).*phi.*((input(t,i)-meanF(:,i))./stdF(:,i).^2);
%             stdF1(:,i)=stdF1(:,i)+(output1(t)-output(t)).*(yl-output1(t)).*phi.*(input(t,i)-meanF(:,i)).^2./stdF(:,i).^3;
%         end
%         g=[yll',reshape(meanF1',[1,M*n]),reshape(stdF1',[1,M*n])];
%         V=beta1*V+(1-beta1)*g;
%         S=beta2*S+(1-beta2)*(g.*g);
%         V_h=V/(1-beta1^t);
%         S_h=S/(1-beta2^2);
%         gg=alpha*V_h./(sqrt(S)+1e-8);
%         xt=xt-gg;
%         yl=xt(1:M)'
%         meanF=reshape(xt(M+1:(n+1)*M),[n,M])';
%         stdF=reshape(xt((n+1)*M+1:end),[n,M])';
%         t
%         %yl=yl1;meanF=meanF1;stdF=stdF1;
%     end
