function [DA_max_x,DA_max] = DA_comput_max(Z_up,W_low,W_up)
N = length(Z_up);
X = diff(Z_up);
S1 = zeros(1,N-1);
S2 = zeros(1,N-1);
for i = 1:N-1
    S1(i) = sum(W_low([1:i]));
    S2(i) = sum(W_up([N:-1:N+1-i]));
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
k = min(find(D > 10e-5));
% if k = N
%     k = N-1;
% end
DA_max_x = [W_low([1:k]); W_up([k+1:N])];
DA_max = sum(DA_max_x.*Z_up)/sum(DA_max_x);

end