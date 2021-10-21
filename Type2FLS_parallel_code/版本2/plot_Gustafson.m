%SAVEspeedUP_ratio_D = zeros(CORE_NUMBER,DATA_SIZE,RULE_NUMBER); 
%1ά�Ƿ���ĺ������ĸ�����2ά�����ݹ�ģ���ĸ�����3ά�ǹ������ĸ���
%% ����SAVEspeedUP_ratio_D�е����ݻ���Gustafson�����й�ͼ
close all

data=SAVEspeedUP_ratio_D;
RGB=[0,0,0; %��ɫ
    1,0,0;  %��ɫ
    1,0.5,0;
    1,0.5,0.5;
    1,0,0.5;
    0.5,0.2,1;
    0.5,0.5,0.5;];

[K,D,R]=size(SAVEspeedUP_ratio_D) 
%��ú������ĸ��������ݹ�ģ���ĸ������������ĸ����˴�K=24,D=7,R=10
%����ֻѡ�����������ݹ�ģ(1,2,3,4,5,6),������Ϊ2,4,6,...,24(12��),����Ϊ2,4,6,8
threadN=[2:2:24];

Fitting_coefficient = zeros(1,4);

figure('color',[1,1,1])   %���ñ�����ɫΪ��ɫ
% for i = 2:2:8  %������20��40��60��80
%     subplot(2,2,i/2)
for i = 2:2:4  %������20��40��60��80
    subplot(1,2,i/2)
    set(gca,'xtick',[2:2:12]);
    set(gca,'FontSize',20);
    box on %�߽�ӿ�
    hold on
    speedup=[];
    for j=1:D-1
        speedup(j)=data(2*j,j,i);
    end
    
    scatter(threadN(1:D-1),speedup,10,RGB(1,:),'filled');
    xlabel('Thread Nunber')
    ylabel('Speedup Ratio')
    title(['Number of rules:',num2str(i*10)])
       %% �ֱ�������С���˷�������ݵ�
%        A=[ones(1,12);1./threadN]';        
%        Fit_coef=(A'*A)\A'*(1./speedup)';%Fit_coef(1)ָ���з����ı�����Fit_coef(2)ָ���з����ı���
%        Fitting_coefficient(:,j,i/2) = Fit_coef;
    X=linspace(2,2*(D-1),200);
%        Y=zeros(1,length(X));
%        f=@(x)1./(Fit_coef(1)+Fit_coef(2)/x);
    a=(sum((threadN(1:D-1)-ones(1,D-1)).*(threadN(1:D-1)-speedup)))/sum((threadN(1:D-1)-ones(1,D-1)).^2);
    Fitting_coefficient(i/2) = a;
    f=@(n) -(n-1)*a+n;
    for s=1:length(X)
        Y(s)=f(X(s));
    end
    l=plot(X,Y,'Color',RGB(j,:));
    xlim([2,12]);
    %legend([l(1),l(2),l(3),l(4),l(5),l(6)],'d500','d1000','d1500','d2000','d2500','d3000')
    hold off  
end
        