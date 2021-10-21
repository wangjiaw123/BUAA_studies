

year =[2014:2020];    %�������
year_num = length(year);
year_time_interval=[283,365;      %2014�ĵ�283�쵽365��
                    1,365;        %2015
                    1,366;
                    1,365;        %2017
                    1,365;
                    1,365;        %2019
                    1,305;];      %2020

Web_address = 'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.mag.calibrated/data/ss/1sec/';    %��Ҫ����ҳԴ����鿴 �����������              
file_save_path = [pwd,'\Data_save'];    %��ȡ ��ǰ���� �����ļ��еĵ�ַ�������н���'Data_save'�ļ��б�������


for i = 1:year_num
    for j = year_time_interval(i,1):year_time_interval(i,2)
        [month_day,month_] = numDays_to_monthDays(year(i),j);   %���� ��-�� �ַ���������:10��21 -> '1021'.
        datafile = ['mvn_mag_l2_',num2str(year(i)),num2str(j,'%03d'),'ss1s_',num2str(year(i)),month_day,'_v01_r01'];  %���ص��ļ���
        data_path = [Web_address,'/',num2str(year(i)),'/',num2str(month_,'%02d'),'/',datafile,'.sts'] ;  % ��Ӧ�ļ������ص�ַ  
        
        save_path = [file_save_path,'\',num2str(year(i)),'\',num2str(month_,'%02d')];   %����·��
        if ~exist(save_path,'dir')
            mkdir(save_path);
        end
        
        datafile_name = [datafile,'.sts'];
        disp(['***��������,',num2str(year(i)),'�꣬',num2str(month_,'%02d'),'�£��ļ���Ϊ',datafile_name,'***'])
        %weboptions('Timeout',20);
        outfilename=websave([save_path,'\',datafile,'.sts'],data_path);   
        
         if ~isempty(outfilename)
             fprintf('*�ļ�%s���سɹ���\n \n',[datafile,'.sts']);
         else
             fprintf('*�ļ�%s����ʧ�ܣ�\n \n',[datafile,'.sts']);
         end
      
         dos_line=['rename',' ',save_path,'\',datafile,'.sts',' ',datafile,'.txt']  ;    
         system(dos_line)  ;                %�޸��ļ���׺��Ϊ.txt

    end
end

fprintf('�����ļ��Ѿ�������ϣ�\n')