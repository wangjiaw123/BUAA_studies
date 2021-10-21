close all
clear
clc
number = 20;
% rand('state',sum(100*clock));
for i = 10:number
%     nn = i;
%     Z_low = 100*rand(nn,1);
%     Z_up = Z_low+10*rand(nn,1);
%     W_low = 5*rand(nn,1);
%     W_up = Z_low+15*rand(nn,1);
    Z_low1 = load('Z_low');
    Z_up1 = load('Z_up');
    W_low1 = load('W_low');
    W_up1 = load('W_up'); 
    Z_low = Z_low1.Z_low([1:i]);
    Z_up = Z_up1.Z_up([1:i]);
    W_low = W_low1.W_low([1:i]);
    W_up = W_up1.W_up([1:i]);
    
%      t=clock;
%      test_KM
%      time(1,i)=etime(clock,t);
%     
%     t=clock;
%     test_EKM
%     time(2,i)=etime(clock,t);
% 
%     t=clock;
%     test_DA
%     time(3,i)=etime(clock,t);
 
%       t=clock;
%       test_KM_my
%       time(4,i)=etime(clock,t);
      
     t=clock;
     test_PSA
     time(1,i)=etime(clock,t);
    
    t=clock;
    test_approximate_method
    time(2,i)=etime(clock,t);
    

    
%     t=clock;
%     test_iterative_method
%     time(3,i)=etime(clock,t);
%     
%     t=clock;
%     test_change_variable_method
%     time(4,i)=etime(clock,t);
    


end
plot([10:number],time([1,2],10:end))
xlabel('N')
ylabel('time')
% legend('KM','EKM','DA','KM_{my}')
% legend('PSA','approximate method','iterative method','change variable method')
legend('PSA','approximate method')