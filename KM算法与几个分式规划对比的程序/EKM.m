%% EKMÀ„∑®
%%
function [EKM_min_x,EKM_min,EKM_max_x,EKM_max] = EKM(Z_low,Z_up,W_low,W_up)
[Z_low,PX] = sort(Z_low);
Z_up = Z_up(PX);
W_low = W_low(PX);
W_up = W_up(PX);
[EKM_min_x,EKM_min] = EKM_comput_min(Z_low,W_low,W_up);
[EKM_max_x,EKM_max] = EKM_comput_max(Z_up,W_low,W_up);
end


