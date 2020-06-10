figure('units','normalized','outerposition',[0 0 1 1])
plots = floor(linspace(1,Nt,10));
for i = 1:size(plots,2)
    % ---------------------------------------------------------------------
    subplot(4,size(plots,2),i);
    subimage(M(:,:,plots(i)))
    %pcolor(rot90(M(:,:,plots(i))')); shading interp;
    %axis square; caxis([0 1]);
    % ---------------------------------------------------------------------
    subplot(4,size(plots,2),i+size(plots,2));
    subimage(M_tilde(:,:,plots(i)))
    %pcolor(rot90(M_tilde(:,:,plots(i))')); shading interp;
    %axis square; caxis([0 1]);
    % ---------------------------------------------------------------------
    subplot(4,size(plots,2),i+2*size(plots,2))
    subimage(M_tilde_moving(:,:,plots(i)))
    %pcolor(rot90(M_tilde_moving(:,:,plots(i))')); shading interp;
    %axis square; caxis([0 1]);
    % ---------------------------------------------------------------------
    subplot(4,size(plots,2),i+3*size(plots,2)); hold on
%     x = bsxfun(@plus, U_x_pos*V_x_pos, reshape(x_e,[Nx*Ny,1]) );
%     y = bsxfun(@plus, U_y_pos*V_y_pos, reshape(y_e,[Nx*Ny,1]) );
    xx = reshape(x(:,plots(i)),[Nx,Ny]);
    yy = reshape(y(:,plots(i)),[Nx,Ny]);
    xxx=xx(1:7:Nx,1:7:Nx); % 7 is chosen to include both ends of the square for Nx = 50
    yyy=yy(1:7:Nx,1:7:Nx);
    scatter(xxx(:),yyy(:),'.k');
    scatter(xxx(1),yyy(1),'or');
end
axes_subM = zeros(1,4*size(plots,2));
for icount = 1:1:4*size(plots,2)
    axes_subM(icount) = subplot(4,size(plots,2),icount); hold on
    rectangle('Position',[0 0 Nx Ny],'EdgeColor','r')
    hold(axes_subM(icount),'on');
    axis(axes_subM(icount),'square');
    axis(axes_subM(icount),'ij');
    box (axes_subM(icount),'on');
    set (axes_subM(icount),'XAxisLocation','top');
    xlim([-0.5 1.5]*Nx)
    ylim([-0.5 1.5]*Ny)
    xlabel('x')
    ylabel('y')

end
linkaxes(axes_subM,'xy');
%%
figure
% x = bsxfun(@plus, U_x_pos*V_x_pos, reshape(x_e,[Nx*Ny,1]) );
% y = bsxfun(@plus, U_y_pos*V_y_pos, reshape(y_e,[Nx*Ny,1]) );
for i = 1:1:Nt
    xx = reshape(x(:,i),[Nx,Ny]);
    yy = reshape(y(:,i),[Nx,Ny]);
    scatter(xx(:),yy(:),i,'.k'); hold on
    axis square;
end
%% Plot Bases
plotter_UV
%%
figure()
ix = 2;
iy = n_mag;
for i = 1:1:n_mag
    subplot(ix,iy,i);
    imagesc(reshape(    U(:,i),Nx,Nx)); title('U_{w} on Eulerian grid')
    axis square;

    subplot(ix,iy,iy+i);
    imagesc(reshape(U_lag(:,i),Nx,Nx)); title('U_{w} on Morphying grid')
    norm(U_lag(:,i))
    axis square;
end
%%
matrix_subplotter_video