function [ output_arg ,output] = AddOfLargeNumbers( input_args1,input_args2 )
if (length(input_args1) > length(input_args2))
    swap=input_args1;
    input_args1=input_args2;
    input_args2=swap;
end
len1=length(input_args1);
len2=length(input_args2);

output=cell(1,len2+1);
for j=1:len2+1
    output{j}='0';
end
carry_num=0;
for i=1:len2
    if i<=len1
        m=str2num(input_args1(len1-i+1))+str2num(input_args2(len2-i+1))+carry_num;
        output{len2+1-i+1}=num2str(mod(m,10));
        carry_num=(m-mod(m,10))/10;
    elseif ((len1+1)<=i)&&(i<=len2)
        m=str2num(input_args2(len2-i+1))+carry_num;
        output{len2-i+2}=num2str(mod(m,10));
        carry_num=(m-mod(m,10))/10 ;
    end
    if (i==len2)&&(carry_num>=0)
         output{1}=num2str(carry_num);
     end
end

for k=1:len2+1
    if (k==1)&&(output{k}~='0')
        outlen=0;
        break;
    elseif (k~=1)&&(output{k}~='0')
        outlen=k-1;
        break;
    elseif (k==1)&&(output{k}=='0') 
        continue;        
    end
end
outlen;
output_arg=cell(1,len2+1-outlen);
for j1=1:len2+1-outlen
    output_arg{j1}=output{j1+outlen};
end

end

