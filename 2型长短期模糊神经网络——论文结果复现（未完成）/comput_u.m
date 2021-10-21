function [ y ] = comput_u( t )
%u 函数u(t)表示控制输入
if t<250
    y=sin(pi*t/25);
elseif (250<=t) && (t<500)
    y=1;
elseif (500<=t) && (t<750)
    y=-1;
else
    y=0.3*sin(pi*t/25)+0.1*sin(pi*t/32)+0.6*sin(pi*t/10);
end
end

