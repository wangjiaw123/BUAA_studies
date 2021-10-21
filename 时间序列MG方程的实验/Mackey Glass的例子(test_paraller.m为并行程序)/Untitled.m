tic 
spmd
    A=rand(labindex)
end
toc


tic
for i=1:6
    A=rand(i)
end
toc

tic
for i=1:600
          sigma_help=zeros(N,n);
          M1_help=zeros(N,n);
          M2_help=zeros(N,n);
          c1_help=zeros(N,1);
          c2_help=zeros(N,1);
          sigma22=sigma;M22=M2;M12=M1;
          c12=c1;c22=c2;
spmd
 [M1,M2,c1,c2,sigma]=train_sfls_type2(XX_train(labindex,:),YY_train(labindex),M1,M2,sigma,c1,c2,alpha);
end
          for ii=1:core_num
              M1_help=M1_help+cell2mat(M1(ii));
              M2_help=M2_help+cell2mat(M2(ii));
              sigma_help=sigma_help+cell2mat(sigma(ii));
              c1_help=c1_help+cell2mat(c1(ii));
              c2_help=c2_help+cell2mat(c2(ii));
          end
end
toc


M1 = rand(N,n);
M2 = M1 + 0.3;
sigma = rand(N,n);
c1 = rand(N,1);
c2 = c1+0.55;
for i=1:3600/500
    X_train=[X_train;X_train];
    Y_train=[Y_train;Y_train];
end


tic
 [M1,M2,c1,c2,sigma]=train_sfls_type2(X_train,Y_train,M1,M2,sigma,c1,c2,alpha);
toc