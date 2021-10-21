function [ value_min,X_min,value_max,X_max ] = theory9_1( Z_low,Z_up,W_low,W_up )

[ value_min,X_min ] = compute_min(Z_low,W_low,W_up);
[ value_max,X_max ] = compute_max(Z_up,W_low,W_up);

end

