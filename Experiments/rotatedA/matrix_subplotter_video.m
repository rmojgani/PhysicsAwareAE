%%
draw_wait_time = 0.025;
figure_video = figure('units','normalized','outerposition',[0 0 1 1]);
pause(1)
axes_video = axes('Parent',figure_video);
hold(axes_video,'on');
axis(axes_video,'square');
axis(axes_video,'ij');
ylim([-0.5 1.5])
xlim([-0.5 1.5])
xlabel('$x$','Interpreter','latex','FontSize', 19)
ylabel('$y$','Interpreter','latex','FontSize', 19)
% Auxiliary objects 
plot([-2,2]     ,[-2,+2]     ,':k')
plot([-1,2]+0.01,[+2,-1]+0.01,':k')
scatter([0.51,0.51],[0.51,0.51],'*k')
% Creating Video Writer Object
writerObj = VideoWriter([file_name,'MAX',num2str(angle.max),'_n_mag_',num2str(n_mag),'_n_pos_',num2str(n_pos),'.avi']);
% Using the 'Open` method to open the file
open(writerObj);

for i = 1:1:Nt
    [x,y] = Ucoarse_2_xyfine(U_x_pos, V_x_pos, U_y_pos, V_y_pos,...
                                size_x,size_y,size_t,n_pos,...
                                Nx, Ny,...
                                x_e, y_e, x_e_c, y_e_c);

    xx = reshape(x(:,i),[Nx,Ny]);
    yy = reshape(y(:,i),[Nx,Ny]);
    xxx=xx(1:7:end,1:7:end);
    yyy=yy(1:7:end,1:7:end);
    
    scatter(xxx(:)./Nx,yyy(:)./Ny,'.k');
    scatter(xxx(1)./Nx,yyy(1)./Ny,'or');
    
    pause(draw_wait_time)
    drawnow
   
   % Adding the frame to the video object using the 'writeVideo' method
    writeVideo(writerObj,  getframe(gcf) );
    
    scatter(xxx(:)./Nx,yyy(:)./Nx,'.g');

    pause(draw_wait_time)

end
% Closing the file and the object using the 'Close' method
close(writerObj);