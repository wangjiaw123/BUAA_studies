function [ output ] = isLeapYear( year )
%isLeapYear �ж��Ƿ�����
if (mod(year,4)==0 && mod(year,100)~=0)||(mod(year,100)==0 && mod(year,400)==0)
    output = true;
else
    output = false;
end

