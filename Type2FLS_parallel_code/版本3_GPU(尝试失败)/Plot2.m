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
    %title(['规则数为：',num2str(j*20)])
    Fig = figure( 'Units', 'pixels','Name', 'move2', 'NumberTitle', 'off', ...
        'IntegerHandle', 'off');
    AxesH = axes( 'Parent', Fig, 'Xlim', [0 30], 'Ylim', [0 0.5], 'XGrid',...
        'on', 'YGrid', 'on', 'DataAspectRatio', [1 1 1], 'Visible', 'on'); 
    for i=1:6 %数据规模
        %subplot(3,2,i)

        %plot(SAVE_RMSE_serial_F(j,:,i))
        line(AxesH,[1:30],SAVE_RMSE_serial_F(j,:,i),'Color',RGB(1,:),'linewidth',2)
        ylabel('RMSE')
        xlabel('epoch')
        hold on
        for k=2:2:24
            %plot(SAVE_RMSE_parallel_G(k,:,i,j))
            line(AxesH,[1:30],SAVE_RMSE_parallel_G(k,:,i,j),'Color',RGB(k/2+1,:))
        end
        legend('s1','p2','p4','p6','p8','p10','p12','p14','p16'...
                ,'p18','p20','p22','p24')
        hold off
        if (j==1 && i==2)||(j==1 && i==3)
            axis([0 30 0 0.6])
        end
        saveas(Fig,['E:\二型模糊逻辑系统的并行实现\结果\7(HPC3)' Fig.Name],'pdf')
    end

end

























