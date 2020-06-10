%% Contour - FOM vs. Projection
tcount = Nt/4%prediction_ratio*Nt;
fig3 = figure();
digit_precision = 1e1;
colorbar_show='on';% 'on' / 'off'
for varcount = 1:1;%numel(fields)
    top1 = -10;
    bottom1 = 10;
    subplot(3, numel(fields), varcount);
    surf(X, Y, M_Eulerian_XY.(fields{varcount})(:,:,tcount) )
    axis equal, shading flat
%     view([0 90])
    view([22 30])

    title( [(fields{varcount}),', t:',num2str(tcount),'/',num2str(Nt)] )
    
    subplot(3, numel(fields), numel(fields)+varcount);
    surf(X, Y, M_Eulerian_Proj_XY.(fields{varcount})(:,:,tcount) )
    axis equal, shading flat
%     view([0 90])
    view([22 30])

    title( [(fields{varcount}),', t: ',num2str(tcount),'/',num2str(Nt)] )

    subplot(3, numel(fields), 2*numel(fields)+varcount);
    surf(X, Y, M_ALE_Proj_on_stationary_grid_XY.(fields{varcount})(:,:,tcount) )
    axis equal, shading flat
%     view([0 90])
    view([22 30])

    title( [(fields{varcount}),', t: ',num2str(tcount),'/',num2str(Nt)] )
    
        
    bottom    = min( min( min(min(min(M_Eulerian_Proj_XY.(fields{varcount})))),...M_Eulerian_XY
                min(min(min(M_Eulerian_Proj_XY.(fields{varcount})))) ),...
                min(min(min(M_ALE_Proj_on_stationary_grid_XY.(fields{varcount})))) );
    bottom_M.(fields{varcount}) = floor(min(bottom1, bottom)*digit_precision)/digit_precision;
    
    top    = max( max( max(max(max(M_Eulerian_Proj_XY.(fields{varcount})))),...
                max(max(max(M_Eulerian_Proj_XY.(fields{varcount})))) ),...
                max(max(max(M_ALE_Proj_on_stationary_grid_XY.(fields{varcount})))) );
    top_M.(fields{varcount}) = ceil(max(top1, top)*digit_precision)/digit_precision;
    
    lim_colorbar.(fields{varcount}) = [bottom_M.(fields{varcount}) top_M.(fields{varcount})];
        
    subplot(3, numel(fields), varcount);
    colorbar('Visible',colorbar_show,'Limits',lim_colorbar.(fields{varcount}) );
    caxis(lim_colorbar.(fields{varcount}))

    subplot(3, numel(fields), numel(fields)+varcount);
    colorbar('Visible',colorbar_show,'Limits',lim_colorbar.(fields{varcount}) );
    caxis(lim_colorbar.(fields{varcount}))
    
    subplot(3, numel(fields), 2*numel(fields)+varcount);
    colorbar('Visible',colorbar_show,'Limits',lim_colorbar.(fields{varcount}) );
    caxis(lim_colorbar.(fields{varcount}))
end
%%
%%
nnx = 75;
nny = 75;
% HFM
uu = M_Eulerian_XY.rho(:,:,tcount); file_name = 'rho_HFM';

% Projection
% uu = M_Eulerian_Proj_XY.rho(:,:,tcount); file_name = 'rho_Eul_Proj';
% uu = M_ALE_Proj_on_stationary_grid_XY.rho(:,:,tcount); file_name = 'rho_ALE_Proj';

% ROM
% uu = M_Eulerian_ROM_XY.rho(:,:,tcount); file_name = 'rho_Eul_ROM';
% uu = M_ALE_ROM_on_stationary_grid_XY.rho(:,:,tcount); file_name = 'rho_ALE_ROM';



ix = floor(linspace(1,150,nnx));
iy = floor(linspace(1,150,nny));
% [X,Y] = meshgrid( x_fine(ix), t_fine(it) );
% [Yp,Xp] = meshgrid( x_fine(ix), t_fine(it) );
Xp = X(ix,iy);
Yp = Y(ix,iy);
Z = uu(ix,iy);
figure();
surf(Xp,Yp,Z,'LineStyle','none'); hold on; 
xlabel('$x$','fontsize',24,'interpreter','latex')
ylabel('$t$','fontsize',24,'interpreter','latex')
zlabel('$w$','fontsize',24,'interpreter','latex')


XX=Xp'; YY=Yp'; ZZ=Z';
data = [ XX(:) YY(:) ZZ(:) ];
% save file_name data -ASCII
save([file_name,'.dat'],'data','-ascii')
% %% diag
iy = floor(linspace(1,Nt,11));
% z = bsxfun(@times, ones(Nx,1),it)/Nt*tmax;z(:,1)=0;
out =  [ones(Nx,1),linspace(0,1,Nx)',diag(uu)    ];
save([file_name,'_diag.dat'],'out','-ascii')

%% Grid
nnx = 51;
nny = 51;
ix = floor(linspace(1,Nx,nnx));
iy = floor(linspace(1,Ny,nny));

out_gridx =  [ix'*0,linspace(0,1,nnx)',x_morphing_fine(ix,iy,tcount)];
out_gridy =  [iy'*0,linspace(0,1,nny)',y_morphing_fine(ix,iy,tcount)'];
% out_gridx =  [ix'*0,linspace(0,1,nnx)',X(ix,iy)'];
% out_gridy =  [iy'*0,linspace(0,1,nny)',Y(ix,iy)];
save Mgridx out_gridx -ASCII
save Mgridy out_gridy -ASCII

% %%
figure();hold on
plot(out_gridx(:,2),out_gridx(:,3:1:end),'k');
plot(out_gridy(:,3:1:end),out_gridy(:,2),'r');