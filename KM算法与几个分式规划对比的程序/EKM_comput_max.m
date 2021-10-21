%% 
function [EKM_max_x,EKM_max] = EKM_comput_max(Z_up,W_low,W_up)
N = length(Z_up);
k = round(N/1.7);
a = sum(Z_up([1:k]).*W_low([1:k])) + sum(Z_up([k+1:N]).*W_up([k+1:N]));
b = sum(W_low([1:k])) + sum(W_up([k+1:N]));
c = a/b;
kk = max(find(Z_up([1:N-1]) < c));
while kk ~=k
    s = sign(kk - k);
    max_k = max(k,kk);
    min_k = min(k,kk);
    aa = a + s*sum(Z_up([min_k+1:max_k]).*(W_low([min_k+1:max_k])-W_up([min_k+1:max_k])));
    bb = b + s*sum(W_low([min_k+1:max_k])-W_up([min_k+1:max_k]));
    cc = aa/bb;
    c = cc;
    a = aa;
    b = bb;
    k = kk; 
end
EKM_max = c;
EKM_max_x = [W_low([1:k]); W_up([k+1:N])];
end