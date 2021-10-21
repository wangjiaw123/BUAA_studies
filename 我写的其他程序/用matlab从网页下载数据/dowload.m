

year =[2014:2020];    %下载年份
year_num = length(year);
year_time_interval=[283,365;      %2014的第283天到365天
                    1,365;        %2015
                    1,366;
                    1,365;        %2017
                    1,365;
                    1,365;        %2019
                    1,305;];      %2020

Web_address = 'https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/maven.mag.calibrated/data/ss/1sec/';    %需要打开网页源代码查看 这个下载链接              
file_save_path = [pwd,'\Data_save'];    %获取 当前程序 所在文件夹的地址，在其中建立'Data_save'文件夹保存数据


for i = 1:year_num
    for j = year_time_interval(i,1):year_time_interval(i,2)
        [month_day,month_] = numDays_to_monthDays(year(i),j);   %返回 月-日 字符串，例如:10月21 -> '1021'.
        datafile = ['mvn_mag_l2_',num2str(year(i)),num2str(j,'%03d'),'ss1s_',num2str(year(i)),month_day,'_v01_r01'];  %下载的文件名
        data_path = [Web_address,'/',num2str(year(i)),'/',num2str(month_,'%02d'),'/',datafile,'.sts'] ;  % 对应文件的下载地址  
        
        save_path = [file_save_path,'\',num2str(year(i)),'\',num2str(month_,'%02d')];   %保存路径
        if ~exist(save_path,'dir')
            mkdir(save_path);
        end
        
        datafile_name = [datafile,'.sts'];
        disp(['***正在下载,',num2str(year(i)),'年，',num2str(month_,'%02d'),'月，文件名为',datafile_name,'***'])
        %weboptions('Timeout',20);
        outfilename=websave([save_path,'\',datafile,'.sts'],data_path);   
        
         if ~isempty(outfilename)
             fprintf('*文件%s下载成功！\n \n',[datafile,'.sts']);
         else
             fprintf('*文件%s下载失败！\n \n',[datafile,'.sts']);
         end
      
         dos_line=['rename',' ',save_path,'\',datafile,'.sts',' ',datafile,'.txt']  ;    
         system(dos_line)  ;                %修改文件后缀名为.txt

    end
end

fprintf('所有文件已经下载完毕！\n')