%SAVEspeedUP_ratio_D = zeros(CORE_NUMBER,DATA_SIZE,RULE_NUMBER); 
%1ά�Ƿ���ĺ������ĸ�����2ά�����ݹ�ģ���ĸ�����3ά�ǹ������ĸ���
%% ����SAVEspeedUP_ratio_D�е����ݻ���Gustafson�����й�ͼ
close all

data=SAVEspeedUP_ratio_D_MEAN;
RGB=[0,0,0; %��ɫ
    1,0,0;  %��ɫ
    1,0.5,0;
    1,0.5,0.5;
    1,0,0.5;
    0.5,0.2,1;
    0.5,0.5,0.5;];

[K,D,R]=size(SAVEspeedUP_ratio_D_MEAN) 
%��ú������ĸ��������ݹ�ģ���ĸ������������ĸ����˴�K=24,D=7,R=10
%����ֻѡ�����������ݹ�ģ(1,2,3,4,5,6),������Ϊ2,4,6,...,24(12��),����Ϊ2,4,6,8
threadN=[2:2:24];

Fitting_coefficient = zeros(1,4);

figure('color',[1,1,1])   %���ñ�����ɫΪ��ɫ
for i = 2:2:8  %������20��40��60��80
    subplot(2,2,i/2)
    %axis([2,12,0.1 1.5])
    box on %�߽�ӿ�
    hold on
    speedup=[];
    for j=1:D
        speedup(j)=data(2*j,j,i/2);
    end
    
    scatter(threadN(1:D),speedup,10,RGB(1,:),'filled');
    xlabel('Thread Nunber')
    ylabel('Speedup Ratio')
    title(['��������Ϊ��',num2str(i*10)])
       %% �ֱ�������С���˷�������ݵ�
%        A=[ones(1,12);1./threadN]';        
%        Fit_coef=(A'*A)\A'*(1./speedup)';%Fit_coef(1)ָ���з����ı�����Fit_coef(2)ָ���з����ı���
%        Fitting_coefficient(:,j,i/2) = Fit_coef;
         X=linspace(1,2*(D),200);
%        Y=zeros(1,length(X));
%        f=@(x)1./(Fit_coef(1)+Fit_coef(2)/x);
     %a=(sum((threadN(1:D)-ones(1,D)).*(threadN(1:D)-speedup)))/sum((threadN(1:D)-ones(1,D)).^2);
     A=zeros(D,2);S=speedup';
     A(:,1)=[2:2:12];A(:,2)=ones(D,1);
     a=(A'*A)\(A'*S)
    %Fitting_coefficient(i/2) = a(1);
%     f=@(n) -(n-1)*a(1)+n;
    c1=a(1);c2=a(2);
    f=@(n) n*c1+c2;
    for s=1:length(X)
        Y(s)=f(X(s));
    end
    l=plot(X,Y,'Color',RGB(j,:));
    %legend([l(1),l(2),l(3),l(4),l(5),l(6)],'d500','d1000','d1500','d2000','d2500','d3000')
    hold off  
end