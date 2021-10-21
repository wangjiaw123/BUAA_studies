function [ output_arg ] = MultiplicationOfLargeNumbers_AuxiliaryFun1( input1,input2 )
if (length(input1) < length(input2))
    swap=input1;
    input_args1=input2;
    input_args2=swap;
end
len1=length(input1);
len2=length(input2);
output=cell(1,len1+len2);
for j=1:len1+len2
    output{j}='0';
end
input2=str2num(input2);
carry_num=0; 
for i = 1:len1
   if i==len1
       m=(str2num(input1(len1-i+1))*input2+carry_num);
       output{len1+len2-i+1}=num2str(mod(m,10));  
       output{len1+len2-i}=num2str((m-mod(m,10))/10);
   else
        m=(str2num(input1(len1-i+1))*input2+carry_num);
        if length(num2str(m))==1
           carry_num=0 ;  
        elseif length(num2str(m))==2
           carry_num=(m-mod(m,10))/10;
        end
        output{len1+len2-i+1}=num2str(mod(m,10))  ; 
   end
end
outlen=len1+len2;
for k=1:len1+len2
    if (k==1)&&(output{k}~='0')
        outlen=len1+len2;
        break
    elseif (k~=1)&&(output{k}~='0')
        outlen=len1+len2-k+1;
        break
    elseif (k==1)&&(output{k}=='0') 
        continue
    end
end
output_arg=cell(1,outlen);
for j1=1:outlen
    output_arg{j1}=output{j1+len1+len2-outlen};
end

end

