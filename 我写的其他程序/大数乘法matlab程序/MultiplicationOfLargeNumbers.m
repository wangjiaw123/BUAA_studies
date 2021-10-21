function [ output_args ] = MultiplicationOfLargeNumbers( input_args1,input_args2 )
%MultiplicationOfLargeNumbers 大数乘法
%   input_args1,input_args2，output均为字符串
if ischar(input_args1)~=1
    input_args1=num2str(input_args1);
end
if ischar(input_args2)~=1
    input_args2=num2str(input_args2);
end
len1=length(input_args1);
len2=length(input_args2);
SET={'0','1','2','3','4','5','6','7','8','9'};

for i1= 1:len1
    if sum(ismember(SET,input_args1(i1)))==0
        disp('输入的字符串必须全部为数字!')
        return;
    end
end
for i2= 1:len2
    if sum(ismember(SET,input_args2(i2)))==0
        disp('输入的字符串必须全部为数字!')
        return;
    end
end 
out=cell(1,len2);
for k=1:len2
   midu={};
   midu= MultiplicationOfLargeNumbers_AuxiliaryFun1(input_args1,input_args2(len2-k+1)) ;
   len_midu=length(midu);
   if k~=1
      for j=1:k-1
          midu{len_midu+j}='0';
      end
   end
   out{1,k}=midu;  
end

sumout=char(out{1,1})';
if length(out)>=2
    for h=2:length(out)
        hh1=h;
        ss1=char(sumout);
        ss2=char(out{1,h})';
        sumout=AddOfLargeNumbers(char(sumout),char(out{1,h})');
    end
end
 output_args =sumout;

end

