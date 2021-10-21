close all

RGB=[0,0,0; %黑色
    1,0,0;  %红色
    1,0.5,0;
    1,0.5,0.5;
    1,0,0.5;
    0,1,0.5;
    0.5,1,0.5;
    0.5,1,0;
    0,0,1;
    0.5,0.2,1;
    1,1,0;
    1,0,1;
    0.5,0.5,0.5;];

for j=1:4 %规则大小
    %figure(j)
    figure('color',[1,1,1])   %设置背景颜色为白色
    %title(['规则数为：',num2str(j*20)])
    for i=1:6 %数据规模
        subplot(3,2,i)
        %plot(SAVE_RMSE_serial_F(j,:,i))
        plot([1:30],SAVE_RMSE_serial_F(j,:,i),'Color',RGB(1,:),'linewidth',2)
        axis([0 30 0 0.28])
        ylabel('RMSE')
        xlabel('epoch')
        hold on
        for k=2:2:24
            %plot(SAVE_RMSE_parallel_G(k,:,i,j))
            plot([1:30],SAVE_RMSE_parallel_G_MEAN(k,:,i,j),'Color',RGB(k/2+1,:))
            title(['数据规模为：',num2str(i*500)])
        end
        legend('sc1','pc2','pc4','pc6','pc8','pc10','pc12','pc14','pc16'...
                ,'pc18','pc20','pc22','pc24','Location', 'northeastoutside' )
        hold off     
end
%print(figure(j), '-dpng', '-r2000', 'my_figure.png')
end
