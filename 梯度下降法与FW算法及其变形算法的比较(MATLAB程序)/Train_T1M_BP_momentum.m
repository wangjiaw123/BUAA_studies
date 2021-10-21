function [meanF,stdF,yl]=Train_T1M_BP_momentum(input,output,meanF,stdF,yl,gama,eta)

% The function T1TSK is to get parameters of a TSK-M FLS from given inputs,outputs
% input: training inputs
% output: training outputs
% meanF: mean of antecedent mfs F_k^l
% stdF: standard deviation of antecedent mfs F_k^l
% yl: centers of consequent fuzzy sets
% alpha: learning rate

[T,~]=size(input);
[M,n]=size(meanF);
meanF1_v=zeros(M,n);
stdF1_v=zeros(M,n);
yl_v=zeros(M,1);

%% Inference by using Mamdani FLS and tuning parameters by BP 

fl=zeros(M,1); % Firing levels f^l
output1=zeros(T,1);
%% ÌÝ¶ÈÏÂ½µ·¨
    for t=1:T
        % Inference by using TSK FLS
        for l=1:M
            u=1;
            for i=1:n
                u=u*exp(-(input(t,i)-meanF(l,i))^2/(2*stdF(l,i)^2));
            end
            fl(l)=u;
        end
        phi=fl/sum(fl);
        output1(t)=yl'*phi; 
    
        % Tuning parameters by BP
    
        meanF1=zeros(M,n);
        stdF1=zeros(M,n);    

        for i=1:n
            meanF1(:,i)=(output1(t)-output(t)).*(yl-output1(t)).*phi.*((input(t,i)-meanF(:,i))./stdF(:,i).^2);
            stdF1(:,i)=(output1(t)-output(t)).*(yl-output1(t)).*phi.*(input(t,i)-meanF(:,i)).^2./stdF(:,i).^3;
        end
        
        % yl=yl-alpha*(output1(t)-output(t))*phi;
        
        
        
        if t>=2
            meanF1_v=gama*meanF1_v+eta*meanF1;
            stdF1_v=gama*stdF1_v+eta*stdF1; 
            yl_v=gama*yl_v+eta*(output1(t)-output(t))*phi;
            meanF=meanF-meanF1_v;
            stdF=stdF-stdF1_v; 
            yl=yl-yl_v;
        else
            meanF1_v=eta*meanF1;
            stdF1_v=eta*stdF1;
            yl_v=eta*(output1(t)-output(t))*phi;
            meanF=meanF-meanF1_v;
            stdF=stdF-stdF1_v;  
            yl=yl-yl_v;
        end
        %meanF=meanF1;stdF=stdF1;
    end
    
%%      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    