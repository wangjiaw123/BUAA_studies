function output=ST1M(input,meanF,stdF,yl)

% The function T1M is to get outputs of a singleton Mamdani FLS from
% input: given inputs
% meanF: mean of antecedent mfs F_k^l
% stdF: standard deviation of antecedent mfs F_k^l
% yl: centers of consequent fuzzy sets

[T,n]=size(input);
[M,n]=size(meanF);

%% Inference and defuzzification
output=zeros(T,1);
fl=zeros(M,1); % Firing levels f^l
for t=1:T
    for l=1:M
         u=1;
         for i=1:n
             u1=exp(-(input(t,i)-meanF(l,i))^2/(2*stdF(l,i)^2));
             u=u*u1;
         end
         fl(l)=u;
    end
    phi=fl/sum(fl);
    output(t)=yl'*phi;
end


end