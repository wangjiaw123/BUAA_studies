
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
f=a/b;
fx=matlabFunction(f)
fmu1=jacobian(f,mu1);
fmu2=jacobian(f,mu2);
fmu3=jacobian(f,mu3);
fmu4=jacobian(f,mu4);
fmu5=jacobian(f,mu5);
fmu6=jacobian(f,mu6);
fmu7=jacobian(f,mu7);
fmu8=jacobian(f,mu8);
fmu9=jacobian(f,mu9);
fmu10=jacobian(f,mu10);
fmu11=jacobian(f,mu11);
fmu12=jacobian(f,mu12);
fmu13=jacobian(f,mu13);
fmu14=jacobian(f,mu14);
fmu15=jacobian(f,mu15);
fmu16=jacobian(f,mu16);
fmu17=jacobian(f,mu17);
fmu18=jacobian(f,mu18);
fmu19=jacobian(f,mu19);
fmu20=jacobian(f,mu20);
fmu21=jacobian(f,mu21);
fmu22=jacobian(f,mu22);
fmu23=jacobian(f,mu23);
fmu24=jacobian(f,mu24);
fmu25=jacobian(f,mu25);
fmu26=jacobian(f,mu26);
fmu27=jacobian(f,mu27);
fmu28=jacobian(f,mu28);
fmu29=jacobian(f,mu29);
fmu30=jacobian(f,mu30);

fsigma1=jacobian(f,sigma1);
fsigma2=jacobian(f,sigma2);
fsigma3=jacobian(f,sigma3);
fsigma4=jacobian(f,sigma4);
fsigma5=jacobian(f,sigma5);
fsigma6=jacobian(f,sigma6);
fsigma7=jacobian(f,sigma7);
fsigma8=jacobian(f,sigma8);
fsigma9=jacobian(f,sigma9);
fsigma10=jacobian(f,sigma10);
fsigma11=jacobian(f,sigma11);
fsigma12=jacobian(f,sigma12);
fsigma13=jacobian(f,sigma13);
fsigma14=jacobian(f,sigma14);
fsigma15=jacobian(f,sigma15);
fsigma16=jacobian(f,sigma16);
fsigma17=jacobian(f,sigma17);
fsigma18=jacobian(f,sigma18);
fsigma19=jacobian(f,sigma19);
fsigma20=jacobian(f,sigma20);
fsigma21=jacobian(f,sigma21);
fsigma22=jacobian(f,sigma22);
fsigma23=jacobian(f,sigma23);
fsigma24=jacobian(f,sigma24);
fsigma25=jacobian(f,sigma25);
fsigma26=jacobian(f,sigma26);
fsigma27=jacobian(f,sigma27);
fsigma28=jacobian(f,sigma28);
fsigma29=jacobian(f,sigma29);
fsigma30=jacobian(f,sigma30);

fc1=jacobian(f,c1);
fc2=jacobian(f,c2);
fc3=jacobian(f,c3);
fc4=jacobian(f,c4);
fc5=jacobian(f,c5);
fc6=jacobian(f,c6);
fc7=jacobian(f,c7);
fc8=jacobian(f,c8);
fc9=jacobian(f,c9);
fc10=jacobian(f,c10);
Jacob=[fmu1 fmu2 fmu3 fmu4 fmu5 fmu6 fmu7 fmu8 fmu9 fmu10...
    fmu11 fmu12 fmu13 fmu14 fmu15 fmu16 fmu17 fmu18 fmu19 fmu20...
    fmu21 fmu22 fmu23 fmu24 fmu25 fmu26 fmu27 fmu28 fmu29 fmu30...
    fsigma1 fsigma2 fsigma3 fsigma4 fsigma5 fsigma6 fsigma7 fsigma8 fsigma9 fsigma10...
    fsigma11 fsigma12 fsigma13 fsigma14 fsigma15 fsigma16 fsigma17 fsigma18 fsigma19 fsigma20...
    fsigma21 fsigma22 fsigma23 fsigma24 fsigma25 fsigma26 fsigma27 fsigma28 fsigma29 fsigma30...
    fc1 fc2]
fun=matlabFunction(Jacob);












