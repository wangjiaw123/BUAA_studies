function [DA_min_x,DA_min] = DA_comput_min(Z_low,W_low,W_up)
N = length(Z_low);
X = diff(Z_low);
S1 = zeros(1,N-1);
S2 = zeros(1,N-1);
for i = 1:N-1
    S1(i) = sum(W_up([1:i]));
    S2(i) = sum(W_low([N:-1:N+1-i]));
end
Tp = zeros(1,N-1);
Tn = zeros(1,N-1);
for i = 1:N-1
    Tp(i) = sum(X(i).*S1(i));
    Tn(i) = sum(X(N-i).*S2(i));  
end
Dp = zeros(1,N);
Dn = zeros(1,N);
for i = 1:N-1
    Dp(i+1) = sum(Tp([1:i]));
    Dn(i+1) = -sum(Tn([1:N-i])); 
end
D = Dp + Dn;
k = min(find(D > 0));
% if k = N
%     k = N-1;
% end
DA_min_x = [W_up([1:k]); W_low([k+1:N])];
DA_min = sum(DA_min_x.*Z_low)/sum(DA_min_x);
end