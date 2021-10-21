function [ value_min,X_min ] = compute_min( ZZ_low,WW_low,WW_up )
[Zl,wz] = sort(ZZ_low);
W_low=WW_low;
W_up=WW_up;

% for i = 1:length(ZZ_up)
%     W_low(i) = WW_low(wz(i));
%     W_up(i) = WW_up(wz(i));
% end
h = (W_up+W_low)./2;
delta=W_up-h;
Wl=h;
eps = 10e-5;
y1 = sum(Zl'*Wl)/sum(Wl);
k = max(find(Zl(1:length(Zl)-1) <= y1));
for i =1:k
    Wl(i) = h(i) + delta(i);
end
for i =k+1:length(Zl)
    Wl(i) = h(i) - delta(i);
end
y2 = sum(Zl'*Wl)/sum(Wl);

while ~(abs(y1-y2) <= 10e-8)
    y1 = y2;
    k = max(find(Zl(1:length(Zl)-1) <= y1));
    for i =1:k
        Wl(i) = h(i) + delta(i);
    end
    for i =k+1:length(Zl)
        Wl(i) = h(i) - delta(i);
    end
    y2 = sum(Zl'*Wl)/sum(Wl);
end
X_min =Wl;
value_min = y2;


end

