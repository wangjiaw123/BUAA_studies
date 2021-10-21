%SAVEspeedUP_ratio_D = zeros(CORE_NUMBER,DATA_SIZE,RULE_NUMBER); 
%1维是分配的核心数的个数，2维是数据规模数的个数，3维是规则数的个数
%% 利用SAVEspeedUP_ratio_D中的数据画与Amdahl定律有关图
close all

data=SAVEspeedUP_ratio_D_MEAN;
RGB=[0,0,0; %黑色
    1,0,0;  %红色
    1,0.5,0;
    1,0.5,0.5;
    1,0,0.5;
    0.5,0.2,1;
    0.5,0.5,0.5;];

[K,D,R]=size(SAVEspeedUP_ratio_D_MEAN) %24 6 4
%获得核心数的个数、数据规模数的个数、规则数的个数此处K=24,D=7,R=10
%文中只选择了六种数据规模(1,2,3,4,5,6),核心数为2,4,6,...,24(12种),规则为2,4,6,8
threadN=[2:2:24];

Fitting_coefficient = zeros(2,6,4);

figure('color',[1,1,1])   %设置背景颜色为白色
for i = 2:2:8  %规则数20，40，60，80
    subplot(2,2,i/2)
    box on %边界加框
    hold on
    for j=1:D
        speedup=data(2:2:24,j,i/2)';
        scatter(threadN,speedup,10,RGB(j,:),'filled');
        l(j)=plot(threadN,speedup,'color',RGB(j,:))
        
        
        xlabel('Thread Nunber')
        ylabel('Speedup Ratio')
        title(['规则数量为：',num2str(i*10)])
       %% 分别利用最小二乘法拟合数据点
        A=[ones(1,12);1./threadN]';        
        Fit_coef=(A'*A)\A'*(1./speedup)';%Fit_coef(1)指串行分量的比例，Fit_coef(2)指并行分量的比例
        Fitting_coefficient(:,j,i/2) = Fit_coef;
        X=linspace(2,24,200);
        Y=zeros(1,length(X));
        f=@(x)1./(Fit_coef(1)+(1-Fit_coef(1))/x);
        %a=(sum(1./speedup)-sum(1./[2:2:24]))/(12-sum(1./[2:2:24]));
        %a=sum((1./threadN-ones(1,12).*(1./threadN-1./speedup)))/sum((1./threadN-ones(1,12)).^2);
        %Fitting_coefficient(i/2,j) = a;
%         f=@(n)1/(a+(1-a)/n);
        for s=1:length(X)
            Y(s)=f(X(s));
        end
        %l(j)=plot(X,Y,'Color',RGB(j,:));
      end   
       
%     legend('d500','d1000','d1500','d2000','d2500','d3000','Location', 'northeastoutside')
%     for j=1:D
%         speedup=data(2:2:24,j,i/2)';
%         scatter(threadN,speedup,10,RGB(j,:),'filled');
%     end
    legend([l(1),l(2),l(3),l(4),l(5),l(6)],'d500','d1000','d1500','d2000','d2500','d3000','Location', 'northeastoutside')
   
    hold off  
end
        
 