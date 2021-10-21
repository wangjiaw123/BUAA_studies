function [ sheet1,sheet2,cow_location,F,cow_location_save,T,Xvalue] = primal_simplex_algorithm( c1,c0,d1,d0,A,b )
%使用Primal simplex algorithm(主单纯型算法)最大化函数F(x)=(c1*x+c0)/(d1*x+d0),同时满足Ax=b,x>=0. 
%c1,d1都是行向量，x是列向量,A的列数>=行数。F是输出的最大值
%主单纯型算法的实现过程可以看《Fractional Programming》的第62—69页。
[sizeA_r, sizeA_c] = size(A);
if sizeA_r ~= rank(A)
    disp('矩阵A行不满秩，方程出现冗余!')
    return;
end
%% 首先找到矩阵A的极大线性无关组的位置（列位置） 
findbaseA = rref(A);                       %将矩阵化成行最简形
cow_location = [];                         %cow_location是矩阵A的极大线性无关组的列位置
[lenfindbaseA_r,lenfindbaseA_c] = size(findbaseA);
for i = 1:lenfindbaseA_c
    findlocction = find(findbaseA(:,i) == 1);
    if (length(findlocction) == 1) && isempty(cow_location)
       cow_location = [i];
    elseif (length(findlocction) == 1) && (length(cow_location) ~= 0)
        gzB = zeros(lenfindbaseA_r,length(cow_location));
        gzB(findlocction,:) = ones(1,length(cow_location));
        if rank(findbaseA(:,cow_location)-gzB) == rank(findbaseA(:,cow_location))
            cow_location = [cow_location i];%注意cow_location一定是升序排列的
        end    
    end
end

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
%% 初始化主单纯型表sheet1（右半边）
sheet1 = zeros(4+sizeA_r,sizeA_c);
sheet1(1,:) = [1:sizeA_c];          %用于标记列号
sheet1(2:1+sizeA_r,:) = A;          %在sheet1中嵌入矩阵A
sheet1_finalcow=[sizeA_c+1 b' 0 0 0]'; 
sheet1(:,sizeA_c+1)=sheet1_finalcow;%添加这一行与上一行是为了在检查调试时与32行sheet2中的xB进行对比，防止出错，此行可以不要
%% 对嵌入在sheet1中的矩阵A进行行变换(依据cow_location)
num_save=[];
for i = 1:length(cow_location)  
    num_locat = find(sheet1(2:1+sizeA_r,cow_location(i))~= 0)+1;
    if i==1 
       sheet1(num_locat(1),:) = sheet1(num_locat(1),:)./sheet1(num_locat(1),cow_location(i));
       for j = 2:sizeA_r+1   
           if j ~= num_locat(1)
              sheet1(j,:) = sheet1(j,:)-(sheet1(j,cow_location(i))./sheet1(num_locat(1),cow_location(i))).* sheet1(num_locat(1),:);%消元
           end
       end
       num_save=[num_locat(1)];
    end

    if i~=1  
       for k = 1:length(num_locat)
           if ~ismember(num_locat(k),num_save)
              num_save = [num_save num_locat(k)];
              jkk1 = num_locat(k);
              break;
           end
       end
       sheet1(jkk1,:) = sheet1(jkk1,:)./sheet1(jkk1,cow_location(i));
       for j = 2:sizeA_r+1   
           if j ~= jkk1
              sheet1(j,:) = sheet1(j,:)-(sheet1(j,cow_location(i))./sheet1(jkk1,cow_location(i))).* sheet1(jkk1,:);
           end
       end   
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
% %% 如果sheet1最后一行全小于0则结束计算
% if all(sheet1(4+sizeA_r,:) <= 0) 
%    return 
% end
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

%% 进入循环
cow_location_save=[cow_location,enterbase_num,outbase_num];
%~isempty(find(sheet1(sizeA_r+4,:)>0))
T=1;
while ~all(sheet1(sizeA_r+4,:) <= 10e-8)  %设置10e-6的误差，避免正的趋于零导致程序进入无限循环
    [ enterbase_num,outbase_num,sheet1,sheet2,F,cow_location ] = fun1( enterbase_num,sheet1,outbase_num,cow_location,c1,c0,d1,d0,A,b );  
    cow_location_save=[cow_location_save;cow_location,enterbase_num,outbase_num;];
%     T=T+1;
%     if T>21
%         break;
%     end
end
Xvalue = A(:,cow_location)\b;
end

