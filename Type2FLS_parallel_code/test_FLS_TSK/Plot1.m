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
Tmax=40;
    figure('color',[1,1,1])   %设置背景颜色为白色
for j=1:8 %规则大小
    %figure(j)

    %title(['规则数为：',num2str(j*20)])
    subplot(2,4,j)
    %plot(SAVE_RMSE_serial_F(j,:,i))
    plot([1:Tmax],SAVE_RMSE_serial(:,j),'Color',RGB(1,:),'linewidth',2)
    axis([0 Tmax 0 1])
    ylabel('RMSE')
    xlabel('epoch')
    hold on
    for k=1:2:23
        %plot(SAVE_RMSE_parallel_G(k,:,i,j))
        plot([1:Tmax],SAVE_RMSE_parallel(k,:,1,j),'Color',RGB((k+1)/2+1,:))%
        title(['规则数目为：',num2str(j*8)])
    end
    legend('s1','p1','p3','p5','p7','p9','p11','p13','p15'...
              ,'p17','p19','p21','p23','Location', 'northeastoutside' )
    hold off     
end
%print(figure(j), '-dpng', '-r2000', 'my_figure.png')















