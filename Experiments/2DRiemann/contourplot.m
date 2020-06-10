%% Contour - FOM vs. ROM
tcount = 500;%nnt;%prediction_ratio*Nt;
fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
colorbar_show='on';% 'on' / 'off'
for varcount =  1:numel(fields)
    top1 = -10;
    bottom1 = 10;
    
    M2_plot = ( M_Eulerian_ROM_XY.(fields{varcount})(:,:,tcount) - 1*M_Eulerian_XY.(fields{varcount})(:,:,tcount));
    M3_plot = ( M_ALE_ROM_on_stationary_grid_XY.(fields{varcount})(:,:,tcount) - 1*M_Eulerian_XY.(fields{varcount})(:,:,tcount) );

    
    subplot(3, numel(fields), varcount); 
    surf(X, Y, M_Eulerian_XY.(fields{varcount})(:,:,tcount) ); 
colormap(redblue(11))
    %     colormap(flipud(gray))
    axis equal, shading flat
    view([0 90])
    title( [(fields{varcount}),', t: ',num2str(tcount),'/',num2str(Nt)] )
    
    subplot(3, numel(fields), numel(fields)+varcount);
    surf(X, Y, M2_plot )
    axis equal, shading flat
    view([0 90])
    title( [(fields{varcount}),', t: ',num2str(tcount),'/',num2str(Nt)] )

    subplot(3, numel(fields), 2*numel(fields)+varcount);
    surf(X, Y, M3_plot )
    axis equal, shading flat
    view([0 90])
    title( [(fields{varcount}),', t: ',num2str(tcount),'/',num2str(Nt)] )
    
    bottom    = -0.4;
%     min( min( min(min(min(M_Eulerian_XY.(fields{varcount})))),...M_Eulerian_XY
%                 min(min(min(M_Eulerian_XY.(fields{varcount})))) ),...
%                 min(min(min(M_Eulerian_XY.(fields{varcount})))) );
    bottom_M.(fields{varcount}) = floor(min(bottom1, bottom)*digit_precision)/digit_precision;
    
    top    = 0.4;
%     max( max( max(max(max(M_Eulerian_XY.(fields{varcount})))),...
%                 max(max(max(M_Eulerian_XY.(fields{varcount})))) ),...
%                 max(max(max(M_Eulerian_XY.(fields{varcount})))) );
    top_M.(fields{varcount}) =  ceil(max(top1, top)*digit_precision)/digit_precision;
       
    lim_colorbar.(fields{varcount}) = 0.6*[-1 1];%[bottom_M.(fields{varcount}) top_M.(fields{varcount})];
        
    subplot(3, numel(fields), varcount);
    colorbar('Visible',colorbar_show,'Limits',lim_colorbar.(fields{varcount}) );
    caxis(lim_colorbar.(fields{varcount}))

    subplot(3, numel(fields), numel(fields)+varcount);
    colorbar('Visible',colorbar_show,'Limits',lim_colorbar.(fields{varcount}) );
    caxis(lim_colorbar.(fields{varcount}))
%     error_2 = ['$error=$', num2str(norm(M2_plot, 'fro')/Nx/Ny,'%10.1e\n') ];
%     text(0.1,0.1,1,error_2,'Color','white','FontSize',20,'Interpreter','latex')
    error_2 = ['$\| M - {U}{V} \|_{F}=$', num2str(norm(M2_plot, 'fro')/Nx/Ny,'%10.1e\n') ];
    text(0.005,0.1,1,error_2,'Color','Black','FontSize',13,'Interpreter','latex')

    subplot(3, numel(fields), 2*numel(fields)+varcount);
    colorbar('Visible',colorbar_show,'Limits',lim_colorbar.(fields{varcount}) );
    caxis(lim_colorbar.(fields{varcount}))
    
%     error_3 = ['$error=$', num2str(norm(M3_plot, 'fro')/Nx/Ny,'%10.1e\n') ];
%     text(0.1,0.1,1,error_3,'Color','white','FontSize',20,'Interpreter','latex')
    error_3 = ['$\| M - \mathcal{G}^{-1}({U}{V}) \|_{F}=$', num2str(norm(M3_plot, 'fro')/Nx/Ny,'%10.1e\n') ];
    text(0.005,0.1,1,error_3,'Color','Black','FontSize',13,'Interpreter','latex')
end
% saveas(fig3,[num2str(tcount),'.png']);