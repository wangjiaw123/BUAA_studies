clear
close all
clc

load 'averageTSKspeedupRatio.mat';

data=averTSKspeedup;
RGB=[0,0,0; %黑色
    1,0,0;  %红色
    1,0.5,0;
    1,0.5,0.5;
    1,0,0.5;
    0.5,0.2,1;
    0.5,0.5,0.5;
    0.5,0,0;];

threadN=[1:2:23];

Fitting_coefficientA = zeros(2,8);

figure('color',[1,1,1])   %设置背景颜色为白色
subplot(1,2,1)
box on %边界加框
hold on
for i = 1:8  %规则数8,16,24,...,64
    speedup=data(1:2:23,1,i)';
    scatter(threadN,speedup,10,RGB(i,:),'filled');
    xlabel('Thread Nunber')
    ylabel('Speedup Ratio')
    %title(['规则数量为：',num2str(i*10)])
    %% 分别利用最小二乘法拟合数据点
    A=[ones(1,12);1./threadN]';        
    Fit_coef=(A'*A)\A'*(1./speedup)';%Fit_coef(1)指串行分量的比例，Fit_coef(2)指并行分量的比例
    Fitting_coefficientA(:,i) = Fit_coef;
    X=linspace(1,23,200);
    Y=zeros(1,length(X));
    f=@(x)1./(Fit_coef(1)+(Fit_coef(2))/x);
    for s=1:length(X)
        Y(s)=f(X(s));
    end
    l(i)=plot(X,Y,'Color',RGB(i,:))  ;     
 end
 legend([l(1),l(2),l(3),l(4),l(5),l(6),l(7),l(8)],'r8','r16','r24','r32','r40','r48','r56','r64','Location', 'northeastoutside' )
 title('(a)','position',[12.5,-0.9])
 %gtext('(a)不同的规则数量，加速比与线程数的变化');
 hold off  

 
 
subplot(1,2,2)
box on %边界加框
hold on
Fitting_coefficientG = zeros(1,8);
for i = 1:8  %规则数8,16,24,...,64
    speedup=[];
    for j=1:8
        speedup(j)=data(2*j-1,1,j)';
    end
    scatter(threadN(1:8),speedup,10,RGB(1,:),'filled');
    xlabel('Thread Nunber')
    ylabel('Speedup Ratio')
    %title(['规则数量为：',num2str(i*10)])
       %% 分别利用最小二乘法拟合数据点
    X=linspace(1,20,200);
    a=(sum((threadN(1:8)-ones(1,8)).*(threadN(1:8)-speedup)))/sum((threadN(1:8)-ones(1,8)).^2);
    Fitting_coefficient(i) = a;
    f=@(n) -(n-1)*a+n;
    for s=1:length(X)
        Y(s)=f(X(s));
    end
    l=plot(X,Y,'Color',RGB(i,:));
end
title('(b)','position',[10,-0.9])
%gtext('(b)问题规模与线程数等比例增加，加速比与线程数的变化')
hold off  




