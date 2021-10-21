%% DAÀ„∑®
%%
function [DA_min_x,DA_min,DA_max_x,DA_max] = DA(Z_low,Z_up,W_low,W_up)
[Z_low,PX] = sort(Z_low);
Z_up = Z_up(PX);
W_low = W_low(PX);
W_up = W_up(PX);
[DA_min_x,DA_min] = DA_comput_min(Z_low,W_low,W_up);
[DA_max_x,DA_max] = DA_comput_max(Z_up,W_low,W_up);
end








