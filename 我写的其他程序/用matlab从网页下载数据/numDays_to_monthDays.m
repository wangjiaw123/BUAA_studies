function [ output,locat ] = numDays_to_monthDays( year,numDays )
%numDays_to_monthDays 将某明年的某一天转化为 月日 格式,输出为字符串。
%   例如将2020年第275天转换为1001


if isLeapYear(year)
    days_per_month =[0,31,29,31,30,31,30,31,31,30,31,30,31];
else
    days_per_month =[0,31,28,31,30,31,30,31,31,30,31,30,31];
end

locat=0;
for i = 1:12
    if (sum(days_per_month(1:i))<numDays)&&(sum(days_per_month(1:i+1))>=numDays)
        locat = i;
        break
    end
end
rem_ = numDays-sum(days_per_month(1:locat));
output = [num2str(locat,'%02d'),num2str(rem_,'%02d')];    
    
end

