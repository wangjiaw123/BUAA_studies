function [meanF,stdF,yl]=Train_T1M_BP(input,output,meanF,stdF,yl,alpha)

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
        output1(t)=yl'*phi ;
        % Tuning parameters by BP
 
        yl=yl-alpha*(output1(t)-output(t))*phi;
        
        for z=1:M
        for i=1:n
            meanF(z,i)=meanF(z,i)-alpha*(output1(t)-output(t))*(yl(z)-output1(t))*phi(z)*((input(t,i)-meanF(z,i))/stdF(z,i)^2);
            stdF(z,i)=stdF(z,i)-alpha*(output1(t)-output(t))*(yl(z)-output1(t))*phi(z)*(input(t,i)-meanF(z,i))^2/(stdF(z,i)^3);
        end
        end
        
    end
 
    
   
    
    
end
