for i = 5:length(out.hist_z)
    plot3(out.hist_x(1,i-4:i-3), out.hist_x(2,i-4:i-3), out.hist_z(i-4:i-3),'LineWidth',0.5, 'Color',[0 0.4470 0.7410]); xlim([-5 5]); xlim([-5 5]); ylim([-5 5]); zlim([-200 250]); grid on; hold on;
    plot3(out.hist_x(1,i-3:i-2), out.hist_x(2,i-3:i-2), out.hist_z(i-3:i-2),'LineWidth',1, 'Color',[0 0.4470 0.7410]);
    plot3(out.hist_x(1,i-2:i-1), out.hist_x(2,i-2:i-1), out.hist_z(i-2:i-1),'LineWidth',1.5, 'Color',[0 0.4470 0.7410]);
    plot3(out.hist_x(1,i-1:i), out.hist_x(2,i-1:i), out.hist_z(i-1:i),'LineWidth',2, 'Color',[0 0.4470 0.7410]);
    plot3(out.hist_x(1,1:i), out.hist_x(2,1:i), out.hist_z(1:i),'x'); hold off;
    
    view(35,40);
    camproj('perspective');
    drawnow;
    pause(0.05);
end