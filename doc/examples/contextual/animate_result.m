if exist('X0','var')
    X0_len = size(X0,2);
else
    X0_len = 0;
end
for i = 5+X0_len:length(out.hist_z) 
    plot3(out.hist_x(1,i-4:i-3), out.hist_x(2,i-4:i-3), out.hist_z(i-4:i-3),'LineWidth',0.5, 'Color',[0 0.4470 0.7410]);
    xlim([-5 5]); xlim([-5 5]); ylim([-5 5]); zlim([-200 250]); grid on; hold on;
    plot3(out.hist_x(1,i-3:i-2), out.hist_x(2,i-3:i-2), out.hist_z(i-3:i-2),'LineWidth',1, 'Color',[0 0.4470 0.7410]);
    plot3(out.hist_x(1,i-2:i-1), out.hist_x(2,i-2:i-1), out.hist_z(i-2:i-1),'LineWidth',1.5, 'Color',[0 0.4470 0.7410]);
    plot3(out.hist_x(1,i-1:i), out.hist_x(2,i-1:i), out.hist_z(i-1:i),'LineWidth',2, 'Color',[0 0.4470 0.7410]);
    plot3(out.hist_x(1,1:i), out.hist_x(2,1:i), out.hist_z(1:i),'x');
    xlabel('x_1'); ylabel('x_2');
    if X0_len > 0
        plot3(out.hist_x(1,1:X0_len),out.hist_x(2,1:X0_len), out.hist_z(1:X0_len),'o'); 
    end
    hold off;
    view(35,40);
    camproj('perspective');
    drawnow;
    pause(0.01); hold off;
end