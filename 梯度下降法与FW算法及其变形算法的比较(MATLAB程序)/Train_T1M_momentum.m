function [M,sigma,c0,stdF1_v,meanF1_v,yl_v]=Train_T1M_momentum(X,D,...
    M,sigma,c0,gama,eta,stdF1_v,meanF1_v,yl_v)
[L,n]=size(X);
[m,n]=size(M);


meanF1=zeros(m,n);
stdF1=zeros(m,n); 
c_F1=zeros(m,1);
for i=1:L
    U=[];
    for j=1:m
        u=1;
        for t=1:n
            u=u*(gaussmf(X(i,t),[sigma(j,1),M(j,1)]));
        end
        U=[U,u];
    end

    fa=U/sum(U);
    f=fa*c0;
    fa=fa';
    e=D(i)-f;
    
 
    
    for l=1:m
        for k=1:n
            meanF1(l,k)=meanF1(l,k)+e*((X(i,k)-M(l,k))/(sigma(l,k)^2))*(c0(l)-f)*U(l)/sum(U);
            stdF1(l,k)=stdF1(l,k)+e*(((X(i,k)-M(l,k))^2)/(sigma(l,k)^3))*(c0(l)-f)*U(l)/sum(U);
        end
    end
    c_F1=c_F1+e*fa;
    %c0=c0+alpha*e*fa;

end
meanF1_v
meanF1_v=gama*meanF1_v+eta*meanF1;
stdF1_v=gama*stdF1_v+eta*stdF1; 
yl_v=gama*yl_v+eta*c_F1;
M=M+meanF1_v;
sigma=sigma+stdF1_v;
c0=c0+yl_v;
 
end