function [ enterbase_num,outbase_num,sheet1,sheet2,F,cow_location ] = fun1( enterbase_num,sheet1,outbase_num,cow_location,c1,c0,d1,d0,A,b )
%fun1 主单纯型法的辅助程序，进基位置是enterbase_num，出基位置是outbase_num
%sheet1,sheet2都是主单纯型法表格的部分表格
[sizeA_r,sizeA_c] = size(A);
%% 改变矩阵A的极大线性无关组的列位置向量cow_location
a1=find(cow_location == outbase_num); %找到出基的位置
cow_location(a1)=enterbase_num;       %enterbase_num这一列进基
%% 初始化主单纯型表sheet2（左半边）
sheet2 = zeros(1+length(cow_location),3);
sheet2(1:length(cow_location),1) = d1(cow_location)';      %初始化dB（详细说明见书本P67Table3.1.2）
sheet2(1:length(cow_location),2) = c1(cow_location)';      %初始化cB
sheet2(1:length(cow_location),3) = (A(:,cow_location)\b)'; %初始化xB
%sheet2(1:length(cow_location),3) = ((A(:,cow_location)'*A(:,cow_location))*A(:,cow_location)*b)';
sheet2(1+length(cow_location),2) = sheet2(1:length(cow_location),2)'*sheet2(1:length(cow_location),3)+c0;%计算z1
sheet2(1+length(cow_location),1) = sheet2(1:length(cow_location),1)'*sheet2(1:length(cow_location),3)+d0;%计算z2
sheet2(1+length(cow_location),3) = sheet2(1+length(cow_location),2)/sheet2(1+length(cow_location),1);    %计算F
F = sheet2(1+length(cow_location),3);
%% 对嵌入在sheet1中的矩阵A的enterbase_num这一列进行行变换
num_locat = find(sheet1(2:1+sizeA_r,enterbase_num)~= 0)+1;
sheet1(num_locat(1),:) = sheet1(num_locat(1),:)./sheet1(num_locat(1),enterbase_num);
for j = 2:sizeA_r+1   
    if j ~= num_locat(1)
        sheet1(j,:) = sheet1(j,:)-(sheet1(j,enterbase_num)./sheet1(num_locat(1),enterbase_num)).* sheet1(num_locat(1),:);
    end
end
%% 计算sheet1中的倒数第3，2，1行分别表示的xj-zj1,dj-zj2,deataj(=z2*(cj+zj1)-z1*(dj-zj2))
U = A(:,cow_location )\A;
for i = 1:sizeA_c
    Zj1 = sheet2(1:sizeA_r,2)'*U(:,i);
    Zj2 = sheet2(1:sizeA_r,1)'*U(:,i);
    sheet1(2+sizeA_r,i) = c1(i)-Zj1;
    sheet1(3+sizeA_r,i) = d1(i)-Zj2;
    sheet1(4+sizeA_r,i) = sheet2(sizeA_r+1,1)*sheet1(2+sizeA_r,i)-sheet2(sizeA_r+1,2)*sheet1(3+sizeA_r,i);
end
%% 找到哪个进基，哪个出基
enterbase_num1 = find(sheet1(4+sizeA_r,1:sizeA_c) == max(sheet1(4+sizeA_r,setdiff([1:sizeA_c],cow_location,'stable'))));  %找出哪个位置进基
enterbase_num=enterbase_num1(1);  %如果最大值相等的值有多个，则取第一个                     
num_zs1 = find(sheet1(2:1+sizeA_r,enterbase_num) > 0);
num_zs2 = find(min(sheet2(num_zs1,3)./sheet1(1+num_zs1,enterbase_num)) == sheet2(num_zs1,3)./sheet1(1+num_zs1,enterbase_num));   %找出哪个位置出基
num_zs3 = num_zs1(num_zs2);
outbase_num = cow_location(num_zs3(1));

cow_location_save=cow_location;
a1=find(cow_location == outbase_num); %找到出基的位置
cow_location(a1)=enterbase_num; 
if det(sheet1(2:sizeA_r+1,cow_location))==0
    cow_location=cow_location_save;
    return; 
else 
    cow_location=cow_location_save;
end
end
