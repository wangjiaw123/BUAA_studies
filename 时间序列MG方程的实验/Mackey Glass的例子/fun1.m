function [ enterbase_num,outbase_num,sheet1,sheet2,F,cow_location ] = fun1( enterbase_num,sheet1,outbase_num,cow_location,c1,c0,d1,d0,A,b )
%fun1 �������ͷ��ĸ������򣬽���λ����enterbase_num������λ����outbase_num
%sheet1,sheet2�����������ͷ����Ĳ��ֱ��
[sizeA_r,sizeA_c] = size(A);
%% �ı����A�ļ��������޹������λ������cow_location
a1=find(cow_location == outbase_num); %�ҵ�������λ��
cow_location(a1)=enterbase_num;       %enterbase_num��һ�н���
%% ��ʼ���������ͱ�sheet2�����ߣ�
sheet2 = zeros(1+length(cow_location),3);
sheet2(1:length(cow_location),1) = d1(cow_location)';      %��ʼ��dB����ϸ˵�����鱾P67Table3.1.2��
sheet2(1:length(cow_location),2) = c1(cow_location)';      %��ʼ��cB
sheet2(1:length(cow_location),3) = (A(:,cow_location)\b)'; %��ʼ��xB
%sheet2(1:length(cow_location),3) = ((A(:,cow_location)'*A(:,cow_location))*A(:,cow_location)*b)';
sheet2(1+length(cow_location),2) = sheet2(1:length(cow_location),2)'*sheet2(1:length(cow_location),3)+c0;%����z1
sheet2(1+length(cow_location),1) = sheet2(1:length(cow_location),1)'*sheet2(1:length(cow_location),3)+d0;%����z2
sheet2(1+length(cow_location),3) = sheet2(1+length(cow_location),2)/sheet2(1+length(cow_location),1);    %����F
F = sheet2(1+length(cow_location),3);
%% ��Ƕ����sheet1�еľ���A��enterbase_num��һ�н����б任
num_locat = find(sheet1(2:1+sizeA_r,enterbase_num)~= 0)+1;
sheet1(num_locat(1),:) = sheet1(num_locat(1),:)./sheet1(num_locat(1),enterbase_num);
for j = 2:sizeA_r+1   
    if j ~= num_locat(1)
        sheet1(j,:) = sheet1(j,:)-(sheet1(j,enterbase_num)./sheet1(num_locat(1),enterbase_num)).* sheet1(num_locat(1),:);
    end
end
%% ����sheet1�еĵ�����3��2��1�зֱ��ʾ��xj-zj1,dj-zj2,deataj(=z2*(cj+zj1)-z1*(dj-zj2))
U = A(:,cow_location )\A;
for i = 1:sizeA_c
    Zj1 = sheet2(1:sizeA_r,2)'*U(:,i);
    Zj2 = sheet2(1:sizeA_r,1)'*U(:,i);
    sheet1(2+sizeA_r,i) = c1(i)-Zj1;
    sheet1(3+sizeA_r,i) = d1(i)-Zj2;
    sheet1(4+sizeA_r,i) = sheet2(sizeA_r+1,1)*sheet1(2+sizeA_r,i)-sheet2(sizeA_r+1,2)*sheet1(3+sizeA_r,i);
end
%% �ҵ��ĸ��������ĸ�����
enterbase_num1 = find(sheet1(4+sizeA_r,1:sizeA_c) == max(sheet1(4+sizeA_r,setdiff([1:sizeA_c],cow_location,'stable'))));  %�ҳ��ĸ�λ�ý���
enterbase_num=enterbase_num1(1);  %������ֵ��ȵ�ֵ�ж������ȡ��һ��                     
num_zs1 = find(sheet1(2:1+sizeA_r,enterbase_num) > 0);
num_zs2 = find(min(sheet2(num_zs1,3)./sheet1(1+num_zs1,enterbase_num)) == sheet2(num_zs1,3)./sheet1(1+num_zs1,enterbase_num));   %�ҳ��ĸ�λ�ó���
num_zs3 = num_zs1(num_zs2);
outbase_num = cow_location(num_zs3(1));

cow_location_save=cow_location;
a1=find(cow_location == outbase_num); %�ҵ�������λ��
cow_location(a1)=enterbase_num; 
if det(sheet1(2:sizeA_r+1,cow_location))==0
    cow_location=cow_location_save;
    return; 
else 
    cow_location=cow_location_save;
end
end
