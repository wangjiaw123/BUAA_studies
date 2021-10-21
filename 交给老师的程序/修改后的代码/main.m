clear ,clear all
close all,clc
tic
y=rand(3,1);
%% 调用comput_u函数计算u
uu=[];
for t=1:1003
    a=comput_u(t);
    uu=[uu,a];
end
u=uu';
%% 计算y
for t=3:1003
    y(t+1)=(y(t)*y(t-1)*y(t-2)*u(t-1)*(y(t-2)-1)+u(t))/(1+(y(t-2))^2+(y(t-1))^2);
end
%% 画图
y1=y(1:1000);
u1=u(1:1000);
figure
scatter(y1,u1,'filled')
axis([-1.1,1.1,-1.1,1.1]);
xlabel('y')
ylabel('u')
title('u and y')
figure
plot([1:1000],y1)
axis([-1,1005,-5,7]);
xlabel('Time step')
ylabel('real output')
%% 将需要的数据放入矩阵
inputdata=[y(3:1002)';y(2:1001)';y(1:1000)';u(3:1002)';u(2:1001)'];
inputdata=inputdata';
outputdata=[y(4:1003)];
%%  初始化
K=4;
Error=0.000001;
Beta=0.25;
Gama=0.0005;%设置学习率,小于0.001效果好一些
n=length(inputdata(1,:));
% sigma_up=1*rand(n*K,1);
sigma_up=0.3*ones(n*K,1);
Sigma_up= sigma_up;
% sigma_low=rand(n*K,1);
sigma_low=0.1*ones(n*K,1);
Sigma_low=sigma_low;
c=1*rand(n*K,1);
ccc=c;
W=2*rand(4*K,8);
WWW=W;
b=5*rand(4*K,4);
bbb=b;
W5=6*rand(K,2);
WWW5=W5;
Layer4h=rand(K,4);
starLayer4h=Layer4h;
saveCYCLIC_o=rand(K,4);
Tmax=100;
count=0;
input=inputdata(1,:);
output=outputdata(1);
TsaveCYCLIC_f=zeros(4,4,1000);
TsaveCYCLIC_i=zeros(4,4,1000);
TsaveCYCLIC_ct=zeros(4,4,1000);
TsaveCYCLIC_o=zeros(4,4,1000);
TsaveCYCLIC_c=zeros(4,4,1000);
Tsave_Layer4h=zeros(4,4,1000);
Tsave_Layer3F=zeros(4,4,1000);
%% 实现误差反向传播，进行参数更新
 while count<1
       for i=1:length(inputdata(:,1))
%% 计算输入数据的输出
          input=inputdata(i,:);
          output=outputdata(i);
         [ Layer2phi,Layer3F,Layer4h,Layer5y,saveCYCLIC_f,saveCYCLIC_i,saveCYCLIC_c,saveCYCLIC_ct,saveCYCLIC_o,sum1,sum2] = comput_Type2FLNN( input,K,sigma_up,sigma_low,c,W,b,Layer4h,saveCYCLIC_o,W5,Beta );
         for ii1=1:K
             for jj=1:4
                 TsaveCYCLIC_f(ii1,jj,i)=saveCYCLIC_f(ii1,jj); 
                 TsaveCYCLIC_i(ii1,jj,i)=saveCYCLIC_i(ii1,jj);
                 TsaveCYCLIC_ct(ii1,jj,i)=saveCYCLIC_ct(ii1,jj);
                 TsaveCYCLIC_o(ii1,jj,i)=saveCYCLIC_o(ii1,jj);
                 TsaveCYCLIC_c(ii1,jj,i)=saveCYCLIC_c(ii1,jj);
                 Tsave_Layer4h(ii1,jj,i)=Layer4h(ii1,jj);
                 Tsave_Layer3F(ii1,jj,i)=Layer3F(ii1,jj);
             end
         end
         E=0.5*(Layer5y-output)^2;
%% 更新参数W5
         for i1=1:K       
             for j1=1:2
                 if j1==1
                     W5(i1,j1)=W5(i1,j1)-Gama*(Layer5y-output)*(1-Beta)*(Layer4h(i1,1)+Layer4h(i1,2))/sum1;
                 elseif j1==2
                     W5(i1,j1)=W5(i1,j1)-Gama*(Layer5y-output)*Beta*(Layer4h(i1,3)+Layer4h(i1,4))/sum2;
                 end    
             end
         end 
         
       
%% 更新参数b    
      for i2=1:4*K %用到了公式28，35，36，37，38
          if mod(i2,4)==1 && mod(i2,4)==2  %选择参数W5的第一列或第二列
             s1=1;
          else 
             s1=2;
          end
          
          for j2=1:4
              if j2==1
                  if mod(i2,4)==0
                     a1=i2/4;
                  else
                      a1=(i2-mod(i2,4))/4+1;
                  end
                 delta=(Layer5y-output)*W5(a1,s1);
                 a2=0;
                 for k1=1:i
                     delta_f111=delta*TsaveCYCLIC_o(a1,j2,k1)*(1-(compute_tanh(TsaveCYCLIC_c(a1,j2,k1)))^2)*TsaveCYCLIC_c(a1,j2,k1)*TsaveCYCLIC_f(a1,j2,k1)*(1-TsaveCYCLIC_f(a1,j2,k1));
                     a2=a2+delta_f111;  
                 end
                 b(i2,j2)=b(i2,j2)-Gama*a2;
              elseif j2==2
                  if mod(i2,4)==0
                     b1=i2/4;
                  else
                      b1=(i2-mod(i2,4))/4+1;
                  end
                  delta=(Layer5y-output)*W5(b1,s1);                
                  b2=0;
                  for k2=1:i
                      delta_i111=delta*TsaveCYCLIC_o(b1,j2,k2)*(1-(compute_tanh(TsaveCYCLIC_c(b1,j2,k2)))^2)*TsaveCYCLIC_ct(b1,j2,k2)*TsaveCYCLIC_i(b1,j2,k2)*(1-TsaveCYCLIC_i(b1,j2,k2));
                      b2=b2+delta_i111;  
                  end                 
                  b(i2,j2)=b(i2,j2)-Gama*b2;
              elseif j2==3
                  if mod(i2,4)==0
                     c1=i2/4;
                  else
                      c1=(i2-mod(i2,4))/4+1;
                  end
                  delta=(Layer5y-output)*W5(c1,s1);                
                  c2=0;
                  for k3=1:i
                      delta_ct111=delta*TsaveCYCLIC_o(c1,j2,k3)*(1-(compute_tanh(TsaveCYCLIC_c(c1,j2,k3)))^2)*TsaveCYCLIC_i(c1,j2,k3)*(1-(compute_tanh(TsaveCYCLIC_ct(c1,j2,k3)))^2);
                      c2=c2+delta_ct111;
                  end
                  b(i2,j2)=b(i2,j2)-Gama* c2;               
              elseif j2==4
                   if mod(i2,4)==0
                     d1=i2/4;
                  else
                      d1=(i2-mod(i2,4))/4+1;
                  end
                  delta=(Layer5y-output)*W5(d1,s1);                 
                  d2=0;
                  for k4=1:i
                      delta_o111=delta*compute_tanh(TsaveCYCLIC_ct(d1,j2,k4))*TsaveCYCLIC_o(d1,j2,k4)*(1-TsaveCYCLIC_o(d1,j2,k4));
                      d2=d2+delta_o111;
                  end
                  b(i2,j2)=b(i2,j2)-Gama*d2;  
             end
          end
      end
      

%% 更新参数W
      for i3=1:4*K
          if mod(i3,4)==1 && mod(i3,4)==2  %选择参数W5的第一列或第二列
             sl=1;
          else 
             sl=2;
          end
          
          if mod(i3,4)==0
             xc1=4;
          else 
             xc1=mod(i3,4);
          end
          
          for j3=1:8
              if j3==1
                  if mod(i3,4)==0
                     e11=i3/4;
                  else
                      e11=(i3-mod(i3,4))/4+1;
                  end
                 delta=(Layer5y-output)*W5(e11,sl);
%                  delta_f=delta*TsaveCYCLIC_o(e11,j3,i)*(1-(compute_tanh(TsaveCYCLIC_c(e11,j3,i)))^2)*TsaveCYCLIC_c(e11,j3,i)*TsaveCYCLIC_f(e11,j3,i)*(1-TsaveCYCLIC_f(e11,j3,i));
%                  W(i3,j3)=W(i3,j3)-Gama*delta_f*Tsave_Layer3F(e11,(j3+1)/2,i);      
                 delta_f=delta*TsaveCYCLIC_o(e11,xc1,i)*(1-(compute_tanh(TsaveCYCLIC_c(e11,xc1,i)))^2)*TsaveCYCLIC_c(e11,xc1,i)*TsaveCYCLIC_f(e11,xc1,i)*(1-TsaveCYCLIC_f(e11,xc1,i));
                 W(i3,j3)=W(i3,j3)-Gama*delta_f*Tsave_Layer3F(e11,xc1,i);    
              elseif j3==2
                  if mod(i3,4)==0
                     a11=i3/4;
                  else
                     a11=(i3-mod(i3,4))/4+1;
                  end
                 delta=(Layer5y-output)*W5(a11,sl);
                 a22=0;
                 for k1=2:i
%                      delta_f1=delta*TsaveCYCLIC_o(a11,j3/2,k1)*(1-(compute_tanh(TsaveCYCLIC_c(a11,j3/2,k1)))^2)*TsaveCYCLIC_c(a11,j3/2,k1)*TsaveCYCLIC_f(a11,j3/2,k1)*(1-TsaveCYCLIC_f(a11,j3/2,k1));
%                      H=Tsave_Layer4h(a11,j3/2,k1-1);
                     delta_f22=delta*TsaveCYCLIC_o(a11,xc1,k1)*(1-(compute_tanh(TsaveCYCLIC_c(a11,xc1,k1)))^2)*TsaveCYCLIC_c(a11,xc1,k1)*TsaveCYCLIC_f(a11,xc1,k1)*(1-TsaveCYCLIC_f(a11,xc1,k1));
                     H=Tsave_Layer4h(a11,xc1,k1-1);
                     a22=a22+delta_f22*H;
                 end
                 W(i3,j3)=W(i3,j3)-Gama*a22;
              elseif j3==3
                  if mod(i3,4)==0
                     f11=i3/4;
                  else
                      f11=(i3-mod(i3,4))/4+1;
                  end
                 delta=(Layer5y-output)*W5(f11,sl);
%                  delta_i=delta*TsaveCYCLIC_o(f11,j3,i)*(1-(compute_tanh(TsaveCYCLIC_c(f11,j3,i)))^2)*TsaveCYCLIC_ct(f11,j3,i)*TsaveCYCLIC_i(f11,j3,i)*(1-TsaveCYCLIC_i(f11,j3,i));
%                  W(i3,j3)=W(i3,j3)-Gama*delta_i*Tsave_Layer3F(f11,(j3+1)/2,i);   
                 delta_i22=delta*TsaveCYCLIC_o(f11,xc1,i)*(1-(compute_tanh(TsaveCYCLIC_c(f11,xc1,i)))^2)*TsaveCYCLIC_ct(f11,xc1,i)*TsaveCYCLIC_i(f11,xc1,i)*(1-TsaveCYCLIC_i(f11,xc1,i));
                 W(i3,j3)=W(i3,j3)-Gama*delta_i22*Tsave_Layer3F(f11,xc1,i);   
              elseif j3==4
                  if mod(i3,4)==0
                     b11=i3/4;
                  else
                      b11=(i3-mod(i3,4))/4+1;
                  end
                 delta=(Layer5y-output)*W5(b11,sl);
                 b22=0;
                 for k2=2:i
%                      delta_i(i3,j3,k2)=delta*TsaveCYCLIC_o(b11,j3,k2)*(1-(compute_tanh(TsaveCYCLIC_c(b11,j3,k2)))^2)*TsaveCYCLIC_ct(b11,j3,k2)*TsaveCYCLIC_i(b11,j3,k2)*(1-TsaveCYCLIC_i(b11,j3,k2));
%                      H=Tsave_Layer4h(b11,j3/2,k2-1);                    
%                      b22=b22+delta_i(i3,j3,k2)*Tsave_Layer4h(b11,j3/2,k2)*H;  
                     delta_i22=delta*TsaveCYCLIC_o(b11,xc1,k2)*(1-(compute_tanh(TsaveCYCLIC_c(b11,xc1,k2)))^2)*TsaveCYCLIC_ct(b11,xc1,k2)*TsaveCYCLIC_i(b11,xc1,k2)*(1-TsaveCYCLIC_i(b11,xc1,k2));
                     H=Tsave_Layer4h(b11,j3/2,k2-1);                    
                     b22=b22+delta_i22*Tsave_Layer4h(b11,j3/2,k2)*H;
                 end
                 W(i3,j3)=W(i3,j3)-Gama*b22;
              elseif j3==5
                  if mod(i3,4)==0
                     g11=i3/4;
                  else
                      g11=(i3-mod(i3,4))/4+1;
                  end
                 delta=(Layer5y-output)*W5(g11,sl);
%                  delta_ct=delta*TsaveCYCLIC_o(g11,(j3+1)/2,i)*(1-(compute_tanh(TsaveCYCLIC_c(g11,(j3+1)/2,i)))^2)*TsaveCYCLIC_i(g11,(j3+1)/2,i)*(1-(compute_tanh(TsaveCYCLIC_ct(g11,(j3+1)/2,i)))^2);
%                  W(i3,j3)=W(i3,j3)-Gama*delta_ct*Tsave_Layer3F(g11,(j3+1)/2,i); 
                 delta_ct22=delta*TsaveCYCLIC_o(g11,xc1,i)*(1-(compute_tanh(TsaveCYCLIC_c(g11,xc1,i)))^2)*TsaveCYCLIC_i(g11,xc1,i)*(1-(compute_tanh(TsaveCYCLIC_ct(g11,xc1,i)))^2);
                 W(i3,j3)=W(i3,j3)-Gama*delta_ct22*Tsave_Layer3F(g11,xc1,i);  
              elseif j3==6
                  if mod(i3,4)==0
                     c11=i3/4;
                  else
                      c11=(i3-mod(i3,4))/4+1;
                  end
                 delta=(Layer5y-output)*W5(c11,sl);
                 c22=0;
                 for k3=2:i
%                       delta_ct=delta*TsaveCYCLIC_o(c11,j3/2,k3)*(1-(compute_tanh(TsaveCYCLIC_c(c11,j3/2,k3)))^2)*TsaveCYCLIC_i(c11,j3/2,k3)*(1-(compute_tanh(TsaveCYCLIC_ct(c11,j3/2,k3)))^2);                
%                       H=Tsave_Layer4h(c11,j3/2,k3-1);                     
%                       c22=c22+delta_ct*Tsave_Layer4h(c11,j3/2,k3)*H;  
                      delta_ct22=delta*TsaveCYCLIC_o(c11,xc1,k3)*(1-(compute_tanh(TsaveCYCLIC_c(c11,xc1,k3)))^2)*TsaveCYCLIC_i(c11,xc1,k3)*(1-(compute_tanh(TsaveCYCLIC_ct(c11,xc1,k3)))^2);                
                      H=Tsave_Layer4h(c11,xc1,k3-1);                     
                      c22=c22+delta_ct22*Tsave_Layer4h(c11,xc1,k3)*H;  
                 end
                 W(i3,j3)=W(i3,j3)-Gama*c22;
              elseif j3==7
                  if mod(i3,4)==0
                     h11=i3/4;
                  else
                      h11=(i3-mod(i3,4))/4+1;
                  end
                 delta=(Layer5y-output)*W5(h11,sl);
%                  delta_o=delta*compute_tanh(TsaveCYCLIC_ct(h11,(j3+1)/2,i))*TsaveCYCLIC_o(h11,(j3+1)/2,i)*(1-TsaveCYCLIC_o(h11,(j3+1)/2,i));
%                  W(i3,j3)=W(i3,j3)-Gama*delta_ct*Tsave_Layer3F(h11,(j3+1)/2,i);  
                 delta_o22=delta*compute_tanh(TsaveCYCLIC_ct(h11,xc1,i))*TsaveCYCLIC_o(h11,xc1,i)*(1-TsaveCYCLIC_o(h11,xc1,i));
                 W(i3,j3)=W(i3,j3)-Gama*delta_o22*Tsave_Layer3F(h11,xc1,i);  
              elseif j3==8
                 if mod(i3,4)==0
                     d11=i3/4;
                  else
                      d11=(i3-mod(i3,4))/4+1;
                 end
                 delta=(Layer5y-output)*W5(d11,sl);
                 d22=0;
                 for k4=2:i
%                      delta_o(i3,j3,k4)=delta*compute_tanh(TsaveCYCLIC_ct(d11,j3/2,k4))*TsaveCYCLIC_o(d11,j3/2,k4)*(1-TsaveCYCLIC_o(d11,j3/2,k4));
%                      H=Tsave_Layer4h(d11,j3/2,k4-1);                 
%                      d22=d22+delta_o(i3,j3,k4)*Tsave_Layer4h(d11,j3/2,k4)*H; 
                     delta_o22=delta*compute_tanh(TsaveCYCLIC_ct(d11,xc1,k4))*TsaveCYCLIC_o(d11,xc1,k4)*(1-TsaveCYCLIC_o(d11,xc1,k4));
                     H=Tsave_Layer4h(d11,xc1,k4-1);                 
                     d22=d22+delta_o22*Tsave_Layer4h(d11,xc1,k4)*H;  
                 end
                 W(i3,j3)=W(i3,j3)-Gama*d22;
              
              end
          end
      end
      

%% 更新参数sigma_up，参数sigma_low，参数c
      resigma_up=reshape(sigma_up,[n,K])';
      resigma_low=reshape(sigma_low,[n,K])';
      re_c=reshape(c,[n,K])';
      
      for i4=1:K    %用到了公式35，36，37，38，45，46，47，48，49，50，51，52，53，54，55，56，57
           for j4=1:n
              if i4==1
                  ttj=1;
              else
                  ttj=4*(i4-1)+1;
              end
              delta1=(Layer5y-output)*W5(i4,1); 
              delta_f1=delta1*TsaveCYCLIC_o(i4,1,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,1,i)))^2)*TsaveCYCLIC_c(i4,1,i)*TsaveCYCLIC_f(i4,1,i)*(1-TsaveCYCLIC_f(i4,1,i)); 
              delta_i1=delta1*TsaveCYCLIC_o(i4,1,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,1,i)))^2)*TsaveCYCLIC_ct(i4,1,i)*TsaveCYCLIC_i(i4,1,i)*(1-TsaveCYCLIC_i(i4,1,i));
              delta_ct1=delta1*TsaveCYCLIC_o(i4,1,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,1,i)))^2)*TsaveCYCLIC_i(i4,1,i)*(1-(compute_tanh(TsaveCYCLIC_ct(i4,1,i)))^2);
              delta_o1=delta1*compute_tanh(TsaveCYCLIC_ct(i4,1,i))*TsaveCYCLIC_o(i4,1,i)*(1-TsaveCYCLIC_o(i4,1,i));       
              p1=delta_f1*W(ttj,1)+delta_i1*W(ttj,3)+delta_ct1*W(ttj,5)+delta_o1*W(ttj,7);
              q1=Tsave_Layer3F(i4,1,i)/resigma_up(i4,j4);
              
              delta2=(Layer5y-output)*W5(i4,1); 
              delta_f2=delta2*TsaveCYCLIC_o(i4,2,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,2,i)))^2)*TsaveCYCLIC_c(i4,2,i)*TsaveCYCLIC_f(i4,2,i)*(1-TsaveCYCLIC_f(i4,2,i)); 
              delta_i2=delta2*TsaveCYCLIC_o(i4,2,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,2,i)))^2)*TsaveCYCLIC_ct(i4,2,i)*TsaveCYCLIC_i(i4,2,i)*(1-TsaveCYCLIC_i(i4,2,i));
              delta_ct2=delta2*TsaveCYCLIC_o(i4,2,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,2,i)))^2)*TsaveCYCLIC_i(i4,2,i)*(1-(compute_tanh(TsaveCYCLIC_ct(i4,2,i)))^2);
              delta_o2=delta2*compute_tanh(TsaveCYCLIC_ct(i4,1,i))*TsaveCYCLIC_o(i4,1,i)*(1-TsaveCYCLIC_o(i4,1,i));       
              p2=delta_f2*W(ttj+1,1)+delta_i2*W(ttj+1,3)+delta_ct2*W(ttj+1,5)+delta_o2*W(ttj+1,7);       
              q2=Tsave_Layer3F(i4,2,i)/resigma_up(i4,j4);
              
              delta3=(Layer5y-output)*W5(i4,2); 
              delta_f3=delta3*TsaveCYCLIC_o(i4,3,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,3,i)))^2)*TsaveCYCLIC_c(i4,3,i)*TsaveCYCLIC_f(i4,3,i)*(1-TsaveCYCLIC_f(i4,3,i)); 
              delta_i3=delta3*TsaveCYCLIC_o(i4,3,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,3,i)))^2)*TsaveCYCLIC_ct(i4,3,i)*TsaveCYCLIC_i(i4,3,i)*(1-TsaveCYCLIC_i(i4,3,i));
              delta_ct3=delta3*TsaveCYCLIC_o(i4,3,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,3,i)))^2)*TsaveCYCLIC_i(i4,3,i)*(1-(compute_tanh(TsaveCYCLIC_ct(i4,3,i)))^2);
              delta_o3=delta3*compute_tanh(TsaveCYCLIC_ct(i4,1,i))*TsaveCYCLIC_o(i4,1,i)*(1-TsaveCYCLIC_o(i4,1,i));       
              p3=delta_f3*W(ttj+2,1)+delta_i3*W(ttj+2,3)+delta_ct3*W(ttj+2,5)+delta_o3*W(ttj+2,7);   
              q3=Tsave_Layer3F(i4,3,i)/resigma_up(i4,j4);
              
              delta4=(Layer5y-output)*W5(i4,2); 
              delta_f4=delta4*TsaveCYCLIC_o(i4,4,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,4,i)))^2)*TsaveCYCLIC_c(i4,4,i)*TsaveCYCLIC_f(i4,4,i)*(1-TsaveCYCLIC_f(i4,4,i)); 
              delta_i4=delta4*TsaveCYCLIC_o(i4,4,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,4,i)))^2)*TsaveCYCLIC_ct(i4,4,i)*TsaveCYCLIC_i(i4,4,i)*(1-TsaveCYCLIC_i(i4,4,i));
              delta_ct4=delta4*TsaveCYCLIC_o(i4,4,i)*(1-(compute_tanh(TsaveCYCLIC_c(i4,4,i)))^2)*TsaveCYCLIC_i(i4,4,i)*(1-(compute_tanh(TsaveCYCLIC_ct(i4,4,i)))^2);
              delta_o4=delta4*compute_tanh(TsaveCYCLIC_ct(i4,1,i))*TsaveCYCLIC_o(i4,1,i)*(1-TsaveCYCLIC_o(i4,1,i));       
              p4=delta_f4*W(ttj+3,1)+delta_i4*W(ttj+3,3)+delta_ct4*W(ttj+3,5)+delta_o4*W(ttj+3,7); 
              q4=Tsave_Layer3F(i4,4,i)/resigma_up(i4,j4);
              
              z1=(input(j4)-re_c(i4,j4)).*exp(-(input(j4)-re_c(i4,j4)).^2./(2.*(resigma_up(i4,j4)).^2))./(resigma_up(i4,j4)).^2;
              z2=(input(j4)-re_c(i4,j4)).*exp(-(input(j4)-re_c(i4,j4)).^2./(2.*(resigma_low(i4,j4)).^2))./(resigma_low(i4,j4)).^2;
              z3=(input(j4)-re_c(i4,j4)).*exp(-(input(j4)-re_c(i4,j4)).^2./(2.*(resigma_low(i4,j4)).^2))./(resigma_low(i4,j4)).^2;
              z4=(input(j4)-re_c(i4,j4)).*exp(-(input(j4)-re_c(i4,j4)).^2./(2*(resigma_up(i4,j4)).^2))./(resigma_up(i4,j4)).^2;
              
              re_c(i4,j4)=re_c(i4,j4)-Gama*(p1*q1*z1+p2*q2*z2+p3*q3*z3+p4*q4*z4);
              zz1=(input(j4)-re_c(i4,j4))^2/((resigma_up(i4,j4)).^3)*exp(-(input(j4)-re_c(i4,j4))^2/(2.*(resigma_up(i4,j4)).^2));
              zz2=-(input(j4)-re_c(i4,j4))^2/((resigma_up(i4,j4)).^3)*exp(-(input(j4)-re_c(i4,j4))^2/(2.*(resigma_up(i4,j4)).^2));
              zzz1=(input(j4)-re_c(i4,j4))^2/((resigma_low(i4,j4)).^3)*exp(-(input(j4)-re_c(i4,j4))^2/(2.*(resigma_low(i4,j4)).^2));
              zzz2=-(input(j4)-re_c(i4,j4))^2/((resigma_low(i4,j4)).^3)*exp(-(input(j4)-re_c(i4,j4))^2/(2.*(resigma_low(i4,j4)).^2));
              resigma_up(i4,j4)=resigma_up(i4,j4)-Gama*(p1*q1*zz1+p3*q3*zz2);
              resigma_low(i4,j4)=resigma_low(i4,j4)-Gama*(p2*q2*zzz1+p4*p4*zzz2);
           end
      end
      sigma_up=reshape(resigma_up',[n*K,1]);
      sigma_low=reshape(resigma_up',[n*K,1]);
      c=reshape(re_c',[n*K,1]);  
      

      [ Layer2phi,Layer3F,Layer4h,Layer5y,saveCYCLIC_f,saveCYCLIC_i,saveCYCLIC_c,saveCYCLIC_ct,saveCYCLIC_o,sum1,sum2] = comput_Type2FLNN( input,K,sigma_up,sigma_low,c,W,b,Layer4h,saveCYCLIC_o,W5,Beta );  
      E=(output-Layer5y)^2;
      if E<Error
          flag1=1;
          break;
      else
          flag1=0;
      end 
      
      end
%% 在参数更新之后，再次计算误差 
      if flag1==1;
          break;
      end
      count=count+1     
 end
%% 画出预测图像
ttt=[];
for i=1:1000
    input=inputdata(i,:);
    output=outputdata(i);
    [ Layer2phi,Layer3F,Layer4h,Layer5y,saveCYCLIC_f,saveCYCLIC_i,saveCYCLIC_c,saveCYCLIC_ct,saveCYCLIC_o,sum1,sum2] = comput_Type2FLNN( input,K,sigma_up,sigma_low,c,W,b,Layer4h,saveCYCLIC_o,W5,Beta );
    ttt=[ttt Layer5y];
end
hold on
plot(ttt)
legend('原始输出','预测输出')
toc
 






