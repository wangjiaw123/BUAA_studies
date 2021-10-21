function [ sheet1,sheet2,cow_location,F,cow_location_save,T,Xvalue] = primal_simplex_algorithm( c1,c0,d1,d0,A,b )
%ʹ��Primal simplex algorithm(���������㷨)��󻯺���F(x)=(c1*x+c0)/(d1*x+d0),ͬʱ����Ax=b,x>=0. 
%c1,d1������������x��������,A������>=������F����������ֵ
%���������㷨��ʵ�ֹ��̿��Կ���Fractional Programming���ĵ�62��69ҳ��
[sizeA_r, sizeA_c] = size(A);
if sizeA_r ~= rank(A)
    disp('����A�в����ȣ����̳�������!')
    return;
end
%% �����ҵ�����A�ļ��������޹����λ�ã���λ�ã� 
findbaseA = rref(A);                       %�����󻯳��������
cow_location = [];                         %cow_location�Ǿ���A�ļ��������޹������λ��
[lenfindbaseA_r,lenfindbaseA_c] = size(findbaseA);
for i = 1:lenfindbaseA_c
    findlocction = find(findbaseA(:,i) == 1);
    if (length(findlocction) == 1) && isempty(cow_location)
       cow_location = [i];
    elseif (length(findlocction) == 1) && (length(cow_location) ~= 0)
        gzB = zeros(lenfindbaseA_r,length(cow_location));
        gzB(findlocction,:) = ones(1,length(cow_location));
        if rank(findbaseA(:,cow_location)-gzB) == rank(findbaseA(:,cow_location))
            cow_location = [cow_location i];%ע��cow_locationһ�����������е�
        end    
    end
end

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
%% ��ʼ���������ͱ�sheet1���Ұ�ߣ�
sheet1 = zeros(4+sizeA_r,sizeA_c);
sheet1(1,:) = [1:sizeA_c];          %���ڱ���к�
sheet1(2:1+sizeA_r,:) = A;          %��sheet1��Ƕ�����A
sheet1_finalcow=[sizeA_c+1 b' 0 0 0]'; 
sheet1(:,sizeA_c+1)=sheet1_finalcow;%�����һ������һ����Ϊ���ڼ�����ʱ��32��sheet2�е�xB���жԱȣ���ֹ�������п��Բ�Ҫ
%% ��Ƕ����sheet1�еľ���A�����б任(����cow_location)
num_save=[];
for i = 1:length(cow_location)  
    num_locat = find(sheet1(2:1+sizeA_r,cow_location(i))~= 0)+1;
    if i==1 
       sheet1(num_locat(1),:) = sheet1(num_locat(1),:)./sheet1(num_locat(1),cow_location(i));
       for j = 2:sizeA_r+1   
           if j ~= num_locat(1)
              sheet1(j,:) = sheet1(j,:)-(sheet1(j,cow_location(i))./sheet1(num_locat(1),cow_location(i))).* sheet1(num_locat(1),:);%��Ԫ
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
%% ����sheet1�еĵ�����3��2��1�зֱ��ʾ��xj-zj1,dj-zj2,deataj(=z2*(cj+zj1)-z1*(dj-zj2))
U = A(:,cow_location )\A;
for i = 1:sizeA_c
    Zj1 = sheet2(1:sizeA_r,2)'*U(:,i);
    Zj2 = sheet2(1:sizeA_r,1)'*U(:,i);
    sheet1(2+sizeA_r,i) = c1(i)-Zj1;
    sheet1(3+sizeA_r,i) = d1(i)-Zj2;
    sheet1(4+sizeA_r,i) = sheet2(sizeA_r+1,1)*sheet1(2+sizeA_r,i)-sheet2(sizeA_r+1,2)*sheet1(3+sizeA_r,i);
end
% %% ���sheet1���һ��ȫС��0���������
% if all(sheet1(4+sizeA_r,:) <= 0) 
%    return 
% end
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

%% ����ѭ��
cow_location_save=[cow_location,enterbase_num,outbase_num];
%~isempty(find(sheet1(sizeA_r+4,:)>0))
T=1;
while ~all(sheet1(sizeA_r+4,:) <= 10e-8)  %����10e-6�����������������㵼�³����������ѭ��
    [ enterbase_num,outbase_num,sheet1,sheet2,F,cow_location ] = fun1( enterbase_num,sheet1,outbase_num,cow_location,c1,c0,d1,d0,A,b );  
    cow_location_save=[cow_location_save;cow_location,enterbase_num,outbase_num;];
%     T=T+1;
%     if T>21
%         break;
%     end
end
Xvalue = A(:,cow_location)\b;
end

