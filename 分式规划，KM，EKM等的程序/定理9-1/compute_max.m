function [ value_max,X_max ] = compute_max( ZZ_up,WW_low,WW_up )
[Zl,wz] = sort(ZZ_up);
W_low=WW_low;
W_up=WW_up;

% for i = 1:length(ZZ_up)
%     W_low(i) = WW_low(wz(i));
%     W_up(i) = WW_up(wz(i));
% end
delta = (W_up-W_low)./2;
h=W_up-delta;
Wl=h;

y1 = sum(Zl'*Wl)/sum(Wl);
k = max(find(Zl(1:length(Zl)-1) <= y1));
for i =1:k
    Wl(i) = h(i) - delta(i);
end
for i =k+1:length(Zl)
    Wl(i) = h(i) + delta(i);
end
y2 = sum(Zl'*Wl)/sum(Wl);

while ~(abs(y1-y2) <= 10e-8)
    y1 = y2;
    k = max(find(Zl(1:length(Zl)-1) <= y1));
    for i =1:k
        Wl(i) = h(i) - delta(i);
    end
    for i =k+1:length(Zl)
        Wl(i) = h(i) + delta(i);
    end
    y2 = sum(Zl'*Wl)/sum(Wl);
end
X_max =Wl;
value_max = y2;

end

