%用Frank-Wolf算法调节模糊系统，使其辨识差分方程y(k+1)=0.3*y(k)+0.6*y(k-1)+g[u(k)],
%其中g(u)=0.6sin(pi*u)+0.3sin(3*pi*u)+0.1sin(5*pi*u),u(k)=cos(2*pi*k/25)
%模糊系统是3输入1输出的，使用

clc,clear
close
clear all
tic
%% 1 产生数据
y=rand(2,1);
u(1)=sin(2*pi/200);
g(1)=0.6*sin(pi*u(1))+0.3*sin(3*pi*u(1))+0.1*sin(5*pi*u(1));
for k=2:402
    u(k)=sin(2*pi*k/200);
    g(k)=0.6*sin(pi*u(k))+0.3*sin(3*pi*u(k))+0.1*sin(5*pi*u(k));
    y(k+1)=0.3*y(k)+0.6*y(k-1)+g(k);
end
plot(y)

M=6;   %设置模糊系统规则个数
xtrain=[y(2:201) y(1:200) g(2:201)'];   %训练模糊系统参数的数据
%[m,n]=size(xtrain);
ytrain=y(3:202);

%% 2 利用MATLAB符号计算来求误差函数（平方误差函数）的偏导数
syms x1 x2 x3 mu1 mu2 mu3 mu4 mu5 mu6 mu7 mu8 mu9 mu10 mu11 mu12 mu13 mu14 ...
    mu15 mu16 mu17 mu18 mu19 mu20 mu21 mu22 mu23 mu24 mu25 mu26 mu27 mu28 ...
    mu29 mu30 sigma1 sigma2 sigma3 sigma4 sigma5 sigma6 sigma7 sigma8 sigma9...
    sigma10 sigma11 sigma12 sigma13 sigma14 sigma15 sigma16 sigma17 sigma18 ...
    sigma19 sigma20 sigma21 sigma22 sigma23 sigma24 sigma25 sigma26 sigma27...
    sigma28 sigma29 sigma30 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10
a= c1*exp(-((x1-mu1)/sigma1)^2)*exp(-((x2-mu2)/sigma2)^2)*exp(-((x3-mu3)/sigma3)^2)...
  +c2*exp(-((x1-mu4)/sigma4)^2)*exp(-((x2-mu5)/sigma5)^2)*exp(-((x3-mu6)/sigma6)^2)...
  +c3*exp(-((x1-mu7)/sigma7)^2)*exp(-((x2-mu8)/sigma8)^2)*exp(-((x3-mu9)/sigma9)^2)...
  +c4*exp(-((x1-mu10)/sigma10)^2)*exp(-((x2-mu11)/sigma11)^2)*exp(-((x3-mu12)/sigma12)^2)...
  +c5*exp(-((x1-mu13)/sigma13)^2)*exp(-((x2-mu14)/sigma14)^2)*exp(-((x3-mu15)/sigma15)^2)...
  +c6*exp(-((x1-mu16)/sigma16)^2)*exp(-((x2-mu17)/sigma17)^2)*exp(-((x3-mu18)/sigma18)^2)...
  +c7*exp(-((x1-mu19)/sigma19)^2)*exp(-((x2-mu20)/sigma20)^2)*exp(-((x3-mu21)/sigma21)^2)...
  +c8*exp(-((x1-mu22)/sigma22)^2)*exp(-((x2-mu23)/sigma23)^2)*exp(-((x3-mu24)/sigma24)^2)...
  +c9*exp(-((x1-mu25)/sigma25)^2)*exp(-((x2-mu26)/sigma26)^2)*exp(-((x3-mu27)/sigma27)^2)...
  +c10*exp(-((x1-mu28)/sigma28)^2)*exp(-((x2-mu29)/sigma29)^2)*exp(-((x3-mu30)/sigma30)^2);

b= exp(-((x1-mu1)/sigma1)^2)*exp(-((x2-mu2)/sigma2)^2)*exp(-((x3-mu3)/sigma3)^2)...
  +exp(-((x1-mu4)/sigma4)^2)*exp(-((x2-mu5)/sigma5)^2)*exp(-((x3-mu6)/sigma6)^2)...
  +exp(-((x1-mu7)/sigma7)^2)*exp(-((x2-mu8)/sigma8)^2)*exp(-((x3-mu9)/sigma9)^2)...
  +exp(-((x1-mu10)/sigma10)^2)*exp(-((x2-mu11)/sigma11)^2)*exp(-((x3-mu12)/sigma12)^2)...
  +exp(-((x1-mu13)/sigma13)^2)*exp(-((x2-mu14)/sigma14)^2)*exp(-((x3-mu15)/sigma15)^2)...
  +exp(-((x1-mu16)/sigma16)^2)*exp(-((x2-mu17)/sigma17)^2)*exp(-((x3-mu18)/sigma18)^2)...
  +exp(-((x1-mu19)/sigma19)^2)*exp(-((x2-mu20)/sigma20)^2)*exp(-((x3-mu21)/sigma21)^2)...
  +exp(-((x1-mu22)/sigma22)^2)*exp(-((x2-mu23)/sigma23)^2)*exp(-((x3-mu24)/sigma24)^2)...
  +exp(-((x1-mu25)/sigma25)^2)*exp(-((x2-mu26)/sigma26)^2)*exp(-((x3-mu27)/sigma27)^2)...
  +exp(-((x1-mu28)/sigma28)^2)*exp(-((x2-mu29)/sigma29)^2)*exp(-((x3-mu30)/sigma30)^2);
F=a/b;
ff=matlabFunction(F);
f=0;
for k=10:11  %%%%%%%%%%%注意此处仅用了*对训练数据
    aaafff=subs(ff,{x1,x2,x3},{xtrain(k,1),xtrain(k,2),xtrain(k,3)});
    xhf=(ytrain(k)-aaafff)^2;%计算平方误差的符号函数
    f=f+xhf;
end
f=matlabFunction(f);%将符号函数转化为函数句柄

f_text=char(f);
fid1=fopen('test1.txt','w');
fprintf(fid1,'%s/n',f_text);
fclose(fid1);

Jacob=jacobian(f,[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,mu1,mu2,mu3,mu4,mu5,mu6,...
    mu7,mu8,mu9,mu10,mu11,mu12,mu13,mu14,mu15,mu16,mu17,mu18,mu19,mu20,mu21,...
    mu22,mu23,mu24,mu25,mu26,mu27,mu28,mu29,mu30,sigma1,sigma2,sigma3,sigma4,...
    sigma5,sigma6,sigma7,sigma8,sigma9,sigma10,sigma11,sigma12,sigma13,sigma14,...
    sigma15,sigma16,sigma17,sigma18,sigma19,sigma20,sigma21,sigma22,sigma23,...
    sigma24,sigma25,sigma26,sigma27,sigma28,sigma29,sigma30]);%求得偏导数矩阵
Jacob_text=char(Jacob)
fid2=fopen('test2.txt','w');
fprintf(fid2,'%s/n',Jacob_text);
fclose(fid2);

funJacob=matlabFunction(Jacob)

%% 3 Frank-Wolf算法
epsilon=10;
x0=10*rand(1,70);
x0_save=x0;
Tmax=1000;
df=funJacob(x0(1),x0(2),x0(3),x0(4),x0(5),x0(6),x0(7),x0(8),x0(9),x0(10),x0(11),x0(12),...
    x0(13),x0(14),x0(15),x0(16),x0(17),x0(18),x0(19),x0(20),x0(21),x0(22),x0(23),x0(24),x0(25),x0(26),...
    x0(27),x0(28),x0(29),x0(30),x0(31),x0(32),x0(33),x0(34),x0(35),x0(36),x0(37),x0(38),x0(39),x0(40),x0(41),...
    x0(42),x0(43),x0(44),x0(45),x0(46),x0(47),x0(48),x0(49),x0(50),x0(51),x0(52),x0(53),x0(54),x0(55),x0(56),...
    x0(57),x0(58),x0(59),x0(60),x0(61),x0(62),x0(63),x0(64),x0(65),x0(66),x0(67),x0(68),x0(69),x0(70));

yk=linprog(df,[],[],[],[],-100*ones(70,1),100*ones(70,1));
%Error=norm(df*(yk-x0));
rho=0.5;
sigma=0.4;
t=0;
fk_save=[];
yk_save=[];
while t<100 %&& Error>epsilon
     %g=feval(gfun,x0); 
     Error=norm(df*(yk-x0'));
     if(Error<epsilon)
         break; 
     end
     m=0; mk=0;
     yk_save=[yk_save yk];
     while(m<20)
         d=yk-x0';
         g=funJacob(x0(1),x0(2),x0(3),x0(4),x0(5),x0(6),x0(7),x0(8),x0(9),x0(10),x0(11),x0(12),...
                       x0(13),x0(14),x0(15),x0(16),x0(17),x0(18),x0(19),x0(20),x0(21),x0(22),x0(23),x0(24),x0(25),x0(26),...
                       x0(27),x0(28),x0(29),x0(30),x0(31),x0(32),x0(33),x0(34),x0(35),x0(36),x0(37),x0(38),x0(39),x0(40),x0(41),...
                       x0(42),x0(43),x0(44),x0(45),x0(46),x0(47),x0(48),x0(49),x0(50),x0(51),x0(52),x0(53),x0(54),x0(55),x0(56),...
                       x0(57),x0(58),x0(59),x0(60),x0(61),x0(62),x0(63),x0(64),x0(65),x0(66),x0(67),x0(68),x0(69),x0(70));
         fvalue0=f(x0(1),x0(2),x0(3),x0(4),x0(5),x0(6),x0(7),x0(8),x0(9),x0(10),x0(11),x0(12),...
                  x0(13),x0(14),x0(15),x0(16),x0(17),x0(18),x0(19),x0(20),x0(21),x0(22),x0(23),x0(24),x0(25),x0(26),...
                  x0(27),x0(28),x0(29),x0(30),x0(31),x0(32),x0(33),x0(34),x0(35),x0(36),x0(37),x0(38),x0(39),x0(40),x0(41),...
                  x0(42),x0(43),x0(44),x0(45),x0(46),x0(47),x0(48),x0(49),x0(50),x0(51),x0(52),x0(53),x0(54),x0(55),x0(56),...
                  x0(57),x0(58),x0(59),x0(60),x0(61),x0(62),x0(63),x0(64),x0(65),x0(66),x0(67),x0(68),x0(69),x0(70));
         x0=x0+rho^m*d';
         fvalue1=f(x0(1),x0(2),x0(3),x0(4),x0(5),x0(6),x0(7),x0(8),x0(9),x0(10),x0(11),x0(12),...
                  x0(13),x0(14),x0(15),x0(16),x0(17),x0(18),x0(19),x0(20),x0(21),x0(22),x0(23),x0(24),x0(25),x0(26),...
                  x0(27),x0(28),x0(29),x0(30),x0(31),x0(32),x0(33),x0(34),x0(35),x0(36),x0(37),x0(38),x0(39),x0(40),x0(41),...
                  x0(42),x0(43),x0(44),x0(45),x0(46),x0(47),x0(48),x0(49),x0(50),x0(51),x0(52),x0(53),x0(54),x0(55),x0(56),...
                  x0(57),x0(58),x0(59),x0(60),x0(61),x0(62),x0(63),x0(64),x0(65),x0(66),x0(67),x0(68),x0(69),x0(70));
          if fvalue1<fvalue0+sigma*rho^m*g*d
             mk=m; 
          break;
          end
          m=m+1;
     end
    x0=x0+rho^mk*d';
    %x0=x0';
    ccc=x0;
    df=funJacob(x0(1),x0(2),x0(3),x0(4),x0(5),x0(6),x0(7),x0(8),x0(9),x0(10),x0(11),x0(12),...
                x0(13),x0(14),x0(15),x0(16),x0(17),x0(18),x0(19),x0(20),x0(21),x0(22),x0(23),x0(24),x0(25),x0(26),...
                x0(27),x0(28),x0(29),x0(30),x0(31),x0(32),x0(33),x0(34),x0(35),x0(36),x0(37),x0(38),x0(39),x0(40),...
                x0(41),x0(42),x0(43),x0(44),x0(45),x0(46),x0(47),x0(48),x0(49),x0(50),x0(51),x0(52),x0(53),x0(54),x0(55),x0(56),...
                x0(57),x0(58),x0(59),x0(60),x0(61),x0(62),x0(63),x0(64),x0(65),x0(66),x0(67),x0(68),x0(69),x0(70));
    yk=linprog(df,[],[],[],[],-100*ones(70,1),100*ones(70,1));
    x0=yk;
    f_yk=f(x0(1),x0(2),x0(3),x0(4),x0(5),x0(6),x0(7),x0(8),x0(9),x0(10),x0(11),x0(12),...
                x0(13),x0(14),x0(15),x0(16),x0(17),x0(18),x0(19),x0(20),x0(21),x0(22),x0(23),x0(24),x0(25),x0(26),...
                x0(27),x0(28),x0(29),x0(30),x0(31),x0(32),x0(33),x0(34),x0(35),x0(36),x0(37),x0(38),x0(39),x0(40),...
                x0(41),x0(42),x0(43),x0(44),x0(45),x0(46),x0(47),x0(48),x0(49),x0(50),x0(51),x0(52),x0(53),x0(54),x0(55),x0(56),...
                x0(57),x0(58),x0(59),x0(60),x0(61),x0(62),x0(63),x0(64),x0(65),x0(66),x0(67),x0(68),x0(69),x0(70));
    fk_save=[fk_save f_yk];
    x0=ccc;
    t=t+1;
end
figure
plot(fk_save)
%% 预测
x0=yk;
double(yk)
f_yc=[];
for k=1:200 
    af1=subs(ff,{c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,mu1,mu2,mu3,mu4,mu5,mu6,...
    mu7,mu8,mu9,mu10,mu11,mu12,mu13,mu14,mu15,mu16,mu17,mu18,mu19,mu20,mu21,...
    mu22,mu23,mu24,mu25,mu26,mu27,mu28,mu29,mu30,sigma1,sigma2,sigma3,sigma4,...
    sigma5,sigma6,sigma7,sigma8,sigma9,sigma10,sigma11,sigma12,sigma13,sigma14,...
    sigma15,sigma16,sigma17,sigma18,sigma19,sigma20,sigma21,sigma22,sigma23,...
    sigma24,sigma25,sigma26,sigma27,sigma28,sigma29,sigma30,x1,x2,x3},...
    {x0(1),x0(2),x0(3),x0(4),x0(5),x0(6),x0(7),x0(8),x0(9),x0(10),x0(11),x0(12),...
     x0(13),x0(14),x0(15),x0(16),x0(17),x0(18),x0(19),x0(20),x0(21),x0(22),x0(23),x0(24),x0(25),x0(26),...
     x0(27),x0(28),x0(29),x0(30),x0(31),x0(32),x0(33),x0(34),x0(35),x0(36),x0(37),x0(38),x0(39),x0(40),...
     x0(41),x0(42),x0(43),x0(44),x0(45),x0(46),x0(47),x0(48),x0(49),x0(50),x0(51),x0(52),x0(53),x0(54),x0(55),x0(56),...
     x0(57),x0(58),x0(59),x0(60),x0(61),x0(62),x0(63),x0(64),x0(65),x0(66),x0(67),x0(68),x0(69),x0(70),xtrain(k,1),xtrain(k,2),xtrain(k,3)});
     f_yc=[f_yc double(af1)];
end
figure
plot(f_yc)

toc