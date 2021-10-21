function [ Layer2phi,Layer3F,Layer4h,Layer5y,saveCYCLIC_f,saveCYCLIC_i,saveCYCLIC_c,saveCYCLIC_ct,saveCYCLIC_o,sum1,sum2] = comput_Type2FLNN( input,K,sigma_up,sigma_low,c,W,b,h_last,c_last,W5,Beta )
%comput_Type2FLNN �����ݻ��Ĵ������ڻ��Ƶ�����2��ֱ��ģ��������
%Layer2phi�ǵڶ��������
%Layer3F�ǵ����������
%Layer4h�ǵ��Ĳ������
%Layer5y�ǵ���������
%saveCYCLIC_f,saveCYCLIC_i,saveCYCLIC_c,saveCYCLIC_ct,saveCYCLIC_o�������
%       ���Ĳ��������г��ֵģ��ڼ������򴫲�ʱ��Ҫ�õ����������������
%sum1,sum2����������֮���ڼ������򴫲��и���W5ʱ��Ҫ�õ����������������
n=length(input(1,:));
Layer2phi=zeros(n*K,4);

%% ����ڶ������Layer2phi
for i=1:n*K
    rem=mod(i,n);
    if rem==0
        lay1input_x=input(n);
    else
        lay1input_x=input(rem);
    end
    Layer2phi(i,1)=exp(-(lay1input_x-c(i))^2/(2*(sigma_up(i))^2));
    Layer2phi(i,2)=exp(-(lay1input_x-c(i))^2/(2*(sigma_low(i))^2));
    Layer2phi(i,3)=1-exp(-(lay1input_x-c(i))^2/(2*(sigma_low(i))^2));
    Layer2phi(i,4)=1-exp(-(lay1input_x-c(i))^2/(2*(sigma_up(i))^2));
end
%% ������������Layer3F
for i=1:K
     s=n*(i-1);
     for j=1:4
         aa1=1;
         for t=1:n
             aa1=aa1*Layer2phi(t+s,j);
         end 
         Layer3F(i,j)=aa1;
     end
end
%% ������Ĳ����Layer4h
for i=1:K 
    for j=1:4
        net_f=W(4*(i-1)+j,1)*Layer3F(i,j)+W(4*(i-1)+j,2)*h_last(i,j)+b(4*(i-1)+j,1);
        cyclic_f=compute_delta(net_f);
        saveCYCLIC_f(i,j)=cyclic_f;
        net_i=W(4*(i-1)+j,3)*Layer3F(i,j)+W(4*(i-1)+j,4)*h_last(i,j)+b(4*(i-1)+j,2);
        cyclic_i=compute_delta(net_i);
        saveCYCLIC_i(i,j)=cyclic_i;
        net_ct=W(4*(i-1)+j,5)*Layer3F(i,j)+W(4*(i-1)+j,6)*h_last(i,j)+b(4*(i-1)+j,3);
        cyclic_ct=compute_tanh(net_ct);
        saveCYCLIC_ct(i,j)=cyclic_ct;
        net_o=W(4*(i-1)+j,7)*Layer3F(i,j)+W(4*(i-1)+j,8)*h_last(i,j)+b(4*(i-1)+j,4);
        cyclic_o=compute_delta(net_o); 
        saveCYCLIC_o(i,j)=cyclic_o;
        cyclic_c=cyclic_f*c_last(i,j)+cyclic_i*cyclic_ct;
        saveCYCLIC_c(i,j)=cyclic_c;
        Layer4h(i,j)=cyclic_o*compute_tanh(cyclic_c);
    end
end
%% �����������Layer5y
sum1=0;sum2=0;
sum1=sum(Layer4h(:,1)+Layer4h(:,2));
sum2=sum(Layer4h(:,3)+Layer4h(:,4)); 
Layer5y=0;
for i=1:K
    Layer5y=Layer5y+(1-Beta)*(Layer4h(i,1)+Layer4h(i,2))*W5(i,1)/sum1+Beta*(Layer4h(i,3)+Layer4h(i,4))*W5(i,2)/sum2;
end
    
end

