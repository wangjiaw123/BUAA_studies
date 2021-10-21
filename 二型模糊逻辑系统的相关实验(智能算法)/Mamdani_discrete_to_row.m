function [ output,M,n ] = Mamdani_discrete_to_row( M1,M2,sigma,c1,c2 )
[M,n]=size(M1);
output=[reshape(M1,[1,M*n]) reshape(M2,[1,M*n]) reshape(sigma,[1,M*n]) ...
    reshape(c1,[1,M]) reshape(c2,[1,M])];

end






