%% 
function [EKM_min_x,EKM_min] = EKM_comput_min(Z_low,W_low,W_up)
N = length(Z_low);
k = round(N/2.4);
a = sum(Z_low([1:k]).*W_up([1:k])) + sum(Z_low([k+1:N]).*W_low([k+1:N]));
b = sum(W_up([1:k])) + sum(W_low([k+1:N]));
c = a/b;
kk = max(find(Z_low([1:N-1]) < c));
while kk ~=k
    s = sign(kk - k);
    max_k = max(k,kk);
    min_k = min(k,kk);
    aa = a + s*sum(Z_low([min_k+1:max_k]).*(W_up([min_k+1:max_k])-W_low([min_k+1:max_k])));
    bb = b + s*sum(W_up([min_k+1:max_k])-W_low([min_k+1:max_k]));
    cc = aa/bb;
    c = cc;
    a = aa;
    b = bb;
    k = kk; 
end
EKM_min = c;
EKM_min_x = [W_up([1:k]); W_low([k+1:N])];
end