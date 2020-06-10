clear all;

load('FOM2.mat',...
	'gamma','Nx','Ny','plot_skip','x_grid','y_grid','dt','t_max',...
	'WRITE_DIRECTORY');

writeMovie = true;
if(writeMovie)
    vidObj = VideoWriter(strcat(WRITE_DIRECTORY,'movie.avi'),'Uncompressed AVI');
    %vidObj.Quality = 50;
    vidObj.FrameRate = 10;    
    open(vidObj);
end
figure('Color',[1 1 1],'pos',[10   10   1024   768])


START_read = 1;
SKIP_read = 1;
END_read = 100;


for it = START_read:SKIP_read:END_read
    fin = strcat(WRITE_DIRECTORY, sprintf('q.%1d',it));
    
    fid = fopen(fin);
    raw = fread(fid, [Nx,Ny*4],'float32');
    fclose(fid);   
    
	rho = raw(1:Nx,1:Ny);
	
	pcolor(x_grid(2:end-1,2:end-1),y_grid(2:end-1,2:end-1),(rho(2:end-1,2:end-1))); 
	shading flat; axis equal; axis tight;
	caxis([0.13 3.82]);	
    drawnow

    if(writeMovie)
      currFrame = getframe(gcf);
      writeVideo(vidObj,currFrame);
    end
end

if(writeMovie)
    close(vidObj);
end
