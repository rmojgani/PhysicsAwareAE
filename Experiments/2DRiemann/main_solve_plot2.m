% % % % clear
% % % % load('main_solve_sm_case12.mat')
% % % % % load('main_solve_sm_k=8_Proj.mat')
% % % % load('main_solve_sm_k=2_ROM.mat')
% % % % % load('main_solve_sm_k=12_ROM')
% % % % load('/media/rmojgani/Data/PhD Work/Euler Equation/Eulerian Grid/2DRiemann/saved/config03/main_solve_sm_k=12_ROM.mat')
% % % % load('/media/rmojgani/Data/PhD Work/Euler Equation/Eulerian Grid/2DRiemann/originalnu/main_solve_sm_k=6_ROM.mat')
% % % kM = 2:2:12;
% % % % % load('main_solve.mat')
% % % %%
% % % % load('/media/rmojgani/Data/PhD Work/Euler Equation/Eulerian Grid/2DRiemann/Untitled Folder/main_solve_pred2d5.mat')
% % % % load('/media/rmojgani/Data/PhD Work/Euler Equation/Eulerian Grid/2DRiemann/Untitled Folder/main_solve_solve2d5.mat','M_Eulerian')
% % % lim_colorbar.p= [-0.6000 2.5000];
% % % lim_colorbar.v= [-0.6000 0.9000];
% % % lim_colorbar.u= [-0.6000 0.9000];
% % % lim_colorbar.rho= [-0.1000 2];
%%
% fig1 = figure();
% fields = fieldnames(error.Eul_Proj);
% for filecount = 1:numel(fields)
%     subplot(1,numel(fields),filecount)
%     title(fields{filecount})
%     hold on
%     loglog(kM,  error.Eul_Proj.(fields{filecount}), '.-k', 'DisplayName','$\| U_{w_{Eul}} V_{w_{Eul}} - {w}_{Eul} \|$')
%     loglog(kM,  error.Eul_ROM.(fields{filecount}), 'o--k', 'DisplayName','$\| U_{w_{Eul}} a_{w_{Eul}} - {w}_{Eul} \|$')
%     loglog(kM,  error.ALE_Proj_PALE_2_Eul.(fields{filecount}), 'o-m', 'DisplayName','$\| P(U_{w_{ALE}}) V_{w_{ALE}} - {w}_{Eul} \|$')
%     loglog(kM,  error.ALE_ROM_PALE_2_Eul.(fields{filecount}), 'o--c', 'DisplayName','$\| P(U_{w_{ALE}}) a_{w_{ALE}} - {w}_{Eul} \|$')
% %     
%     % loglog(kM, error.ALE_Proj_ALE2ALE, '.-r', 'DisplayName','$\| U_{w_{ALE}} V_{w_{ALE}} - {w}_{ALE} \|$')
%     % loglog(kM, error.ALE_Proj_PALE_2_PALE, '*-r', 'DisplayName','$\| P(U_{w_{ALE}}) V_{w_{ALE}} - P({w}_{ALE}) \|$')
%     % loglog(kM, error.ALE_ROM_PALE_2_PALE, '*--g', 'DisplayName','$\| P(U_{w_{ALE}}) a_{w_{ALE}} - P({w}_{ALE}) \|$')
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
%     lgd = legend('Location','northeast');
%     set(lgd,'Interpreter','latex','FontSize',14);
%     xlabel('$k_w$','Interpreter','latex','FontSize',20);
%     ylabel('$error$','Interpreter','latex','FontSize',20);
% %     xlim([0 12])
% %     ylim([1e2 1e4])
%     grid on
% end
% pause(1.005)
%%
fields = fieldnames(M_Eulerian);
nnt = size(M_Eulerian.p,2);
for filecount = 1:numel(fields)
    M_Eulerian_XY.(fields{filecount}) = reshape(  M_Eulerian.(fields{filecount}), [Nx,Ny,nnt]);
    M_Eulerian_Proj_XY.(fields{filecount}) = reshape(  M_Eulerian_Proj.(fields{filecount}), [Nx,Ny,nnt]);
    M_ALE_Proj_on_stationary_grid_XY.(fields{filecount}) = reshape(  M_ALE_Proj_on_stationary_grid.(fields{filecount}), [Nx,Ny,nnt]);
end
% for varcount = 1:numel(fields)
%     M_ALE_Proj_on_stationary_grid_XY.(fields{varcount}) = permute(M_ALE_Proj_on_stationary_grid_XY.(fields{varcount}), [2,1,3] );
% end
for filecount = 1:numel(fields)
    M_Eulerian_XY.(fields{filecount}) = reshape(  M_Eulerian.(fields{filecount}), [Nx,Ny,nnt]);
    M_Eulerian_ROM_XY.(fields{filecount}) = reshape(  M_Eulerian_ROM.(fields{filecount}), [Nx,Ny,nnt]);
    M_ALE_ROM_on_stationary_grid_XY.(fields{filecount}) = reshape(  M_ALE_ROM_on_stationary_grid.(fields{filecount}), [Nx,Ny,nnt]);
end
% for varcount = 1:numel(fields)
%     M_ALE_ROM_on_stationary_grid_XY.(fields{varcount}) = permute(M_ALE_ROM_on_stationary_grid_XY.(fields{varcount}), [2,1,3] );
% end
% % % %%
% % % for varcount = 1:numel(fields)
% % %     size_t_plot = 10;
% % %     Ax = zeros(1,size_t_plot);
% % %     tcounter = 1;
% % %     fig_ROM = figure('units','normalized','outerposition',[0 0 1 1]);
% % %     pause(0.5)
% % %     for tcount = floor(linspace(1,prediction_ratio*Nt,size_t_plot))
% % %         Ax(tcounter) = subplot(5,2,tcounter);
% % %         plot(diag(M_Eulerian_XY.(fields{varcount})(:,:,tcount))                   ,'LineWidth',4, 'color',0.75*[1 1 1],'DisplayName', 'FOM'); hold on
% % %         plot(diag(M_ALE_Proj_on_stationary_grid_XY.(fields{varcount})(:,:,tcount)) ,'.:r', 'DisplayName', 'Morphing grid Projection'); hold on
% % %         plot(diag(M_Eulerian_Proj_XY.(fields{varcount})(:,:,tcount))               ,'.-b', 'DisplayName', 'Eulerian grid Projection');
% % %     %     hold off
% % %         title( ['t: ',num2str(tcount),'/',num2str(Nt)] )
% % %         xlabel('Nx')
% % %         ylabel(fields{varcount})
% % % %         ylim(lim_colorbar.(fields{varcount}))
% % %         xlim([0 Nx])
% % %         lgd = legend('Location','best');
% % %         %set(lgd,'Interpreter','latex','FontSize',14);
% % %         tcounter=tcounter+1;
% % %     end
% % %     drawnow
% % %     saveas(fig_ROM,['Proj_line_plot_diag_',fields{varcount},'.png'])
% % %     pause(1)
% % % end
%% Contour - FOM vs. Projection
tcount = nnt;%prediction_ratio*Nt;
fig3 = figure();
digit_precision = 1e1;
colorbar_show='on';% 'on' / 'off'
for varcount = 1:numel(fields)
    top1 = -10;
    bottom1 = 10;
    subplot(3, numel(fields), varcount);
    surf(X, Y, M_Eulerian_XY.(fields{varcount})(:,:,tcount) )
    axis equal, shading flat
    view([0 90])
    title( [(fields{varcount}),', t:',num2str(tcount),'/',num2str(Nt)] )
    
    subplot(3, numel(fields), numel(fields)+varcount);
    surf(X, Y, M_Eulerian_Proj_XY.(fields{varcount})(:,:,tcount) )
    axis equal, shading flat
    view([0 90])
    title( [(fields{varcount}),', t: ',num2str(tcount),'/',num2str(Nt)] )

    subplot(3, numel(fields), 2*numel(fields)+varcount);
    surf(X, Y, M_ALE_Proj_on_stationary_grid_XY.(fields{varcount})(:,:,tcount) )
    axis equal, shading flat
    view([0 90])
    title( [(fields{varcount}),', t: ',num2str(tcount),'/',num2str(Nt)] )
    
        
%     bottom    = min( min( min(min(min(M_Eulerian_Proj_XY.(fields{varcount})))),...M_Eulerian_XY
%                 min(min(min(M_Eulerian_Proj_XY.(fields{varcount})))) ),...
%                 min(min(min(M_ALE_Proj_on_stationary_grid_XY.(fields{varcount})))) );
%     bottom_M.(fields{varcount}) = floor(min(bottom1, bottom)*digit_precision)/digit_precision;
%     
%     top    = max( max( max(max(max(M_Eulerian_Proj_XY.(fields{varcount})))),...
%                 max(max(max(M_Eulerian_Proj_XY.(fields{varcount})))) ),...
%                 max(max(max(M_ALE_Proj_on_stationary_grid_XY.(fields{varcount})))) );
%     top_M.(fields{varcount}) = ceil(max(top1, top)*digit_precision)/digit_precision;
%     
%     lim_colorbar.(fields{varcount}) = [bottom_M.(fields{varcount}) top_M.(fields{varcount})];
%         
%     subplot(3, numel(fields), varcount);
%     colorbar('Visible',colorbar_show,'Limits',lim_colorbar.(fields{varcount}) );
%     caxis(lim_colorbar.(fields{varcount}))
% 
%     subplot(3, numel(fields), numel(fields)+varcount);
%     colorbar('Visible',colorbar_show,'Limits',lim_colorbar.(fields{varcount}) );
%     caxis(lim_colorbar.(fields{varcount}))
%     
%     subplot(3, numel(fields), 2*numel(fields)+varcount);
%     colorbar('Visible',colorbar_show,'Limits',lim_colorbar.(fields{varcount}) );
%     caxis(lim_colorbar.(fields{varcount}))
end
%%
% % % saveas(fig1,[num2str(1),'.png'])
% % % saveas(fig2,[num2str(2),'.png'])
% % % saveas(fig3,[num2str(3),'.png'])
%% Contour - FOM vs. ROM
ifcolormap=[];%'redblue'
tcount = 100;%nnt;%prediction_ratio*Nt;
fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
colorbar_show='on';% 'on' / 'off'
for varcount =  1:numel(fields)
    top1 = -10;
    bottom1 = 10;
    
    M2_plot = abs( M_Eulerian_ROM_XY.(fields{varcount})(:,:,tcount) - 1*M_Eulerian_XY.(fields{varcount})(:,:,tcount));
    M3_plot = abs( M_ALE_ROM_on_stationary_grid_XY.(fields{varcount})(:,:,tcount) - 1*M_Eulerian_XY.(fields{varcount})(:,:,tcount) );

    
    subplot(3, numel(fields), varcount); 
    surf(X, Y, M_Eulerian_XY.(fields{varcount})(:,:,tcount) ); 
    if strcmp(ifcolormap,'redblue')
        colormap(redblue(11))
    else
        colormap(flipud(gray))
    end
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
    
    bottom    = min( min( min(min(min(M_Eulerian_XY.(fields{varcount})))),...M_Eulerian_XY
                min(min(min(M_Eulerian_XY.(fields{varcount})))) ),...
                min(min(min(M_Eulerian_XY.(fields{varcount})))) );
    bottom_M.(fields{varcount}) = floor(min(bottom1, bottom)*digit_precision)/digit_precision;
    
    top    = max( max( max(max(max(M_Eulerian_XY.(fields{varcount})))),...
                max(max(max(M_Eulerian_XY.(fields{varcount})))) ),...
                max(max(max(M_Eulerian_XY.(fields{varcount})))) );
    top_M.(fields{varcount}) =  ceil(max(top1, top)*digit_precision)/digit_precision;
       
    lim_colorbar.(fields{varcount}) = [0 0.6];%[bottom_M.(fields{varcount}) top_M.(fields{varcount})];
        
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

%% Line plot - FOM and ROM
% % % for varcount = 1%:numel(fields)
% % %     size_t_plot = 10;
% % %     Ax = zeros(1,size_t_plot);
% % %     tcounter = 1;
% % %     fig_ROM = figure('units','normalized','outerposition',[0 0 1 1]);
% % %     pause(0.5)
% % %     for tcount = floor(linspace(1,prediction_ratio*Nt,size_t_plot))
% % %         Ax(tcounter) = subplot(5,2,tcounter);
% % %         plot(diag(M_Eulerian_XY.(fields{varcount})(:,:,tcount))                   ,'LineWidth',4, 'color',0.75*[1 1 1],'DisplayName', 'FOM'); hold on
% % %         plot(diag(M_ALE_ROM_on_stationary_grid_XY.(fields{varcount})(:,:,tcount)) ,'.:r', 'DisplayName', 'Morphing grid ROM'); hold on
% % %         plot(diag(M_Eulerian_ROM_XY.(fields{varcount})(:,:,tcount))               ,'.-b', 'DisplayName', 'Eulerian grid ROM');
% % %         hold off
% % %         title( ['t: ',num2str(tcount),'/',num2str(Nt)] )
% % %         xlabel('Nx')
% % %         ylabel(fields{varcount})
% % %         ylim(lim_colorbar.(fields{varcount}))
% % %         xlim([0 Nx])
% % %         lgd = legend('Location','best');
% % %         %set(lgd,'Interpreter','latex','FontSize',14);
% % %         tcounter=tcounter+1;
% % %     end
% % %     drawnow
% % % %     saveas(fig_ROM,['ROM_line_plot_diag_',fields{varcount},'.png'])
% % %     pause(1)
% % % end
%% Line plots - FOM, Projection, ROM
for varcount = 1%:numel(fields)
    size_t_plot = 4;
    Ax = zeros(1,size_t_plot);
    tcounter = 1;
    fig_ROM = figure('units','normalized','outerposition',[0 0 1 1]);
%     pause(0.5)
    for tcount = floor(linspace(1,prediction_ratio*Nt,size_t_plot))
        Ax(tcounter) = subplot(2,2,tcounter);
        hold on
        plot(diag(M_Eulerian_XY.(fields{varcount})(:,:,tcount))                   ,'LineWidth', 4, 'color',0.75*[1 1 1],'DisplayName', 'FOM'); hold on
        
        plot(diag(M_ALE_ROM_on_stationary_grid_XY.(fields{varcount})(:,:,tcount)) ,':g', 'LineWidth', 2, 'DisplayName', 'Morphing grid ROM');
        plot(diag(M_Eulerian_ROM_XY.(fields{varcount})(:,:,tcount))               ,':k', 'LineWidth', 2, 'DisplayName', 'Eulerian grid ROM');

        
        plot(diag(M_ALE_Proj_on_stationary_grid_XY.(fields{varcount})(:,:,tcount)) ,'.-r', 'DisplayName', 'Morphing grid Projection');
        plot(diag(M_Eulerian_Proj_XY.(fields{varcount})(:,:,tcount))               ,'.-b', 'DisplayName', 'Eulerian grid Projection');
        
        title( ['t: ',num2str(tcount),'/',num2str(Nt)] )
        xlabel('Nx')
        ylabel(fields{varcount})
%         ylim(lim_colorbar.(fields{varcount}))
        xlim([0 Nx])
        lgd = legend('Location','best');
        %set(lgd,'Interpreter','latex','FontSize',14);
        tcounter=tcounter+1;
    end
%     drawnow
%     saveas(fig_ROM,['ROM_Proj_line_plot_diag_',fields{varcount},'.png'])
%     pause(1)
end
%% Grid
fig_grid = figure('units','normalized','outerposition',[0 0 1 1]);
size_t_plot = 6; prediction_ratio=1;
Ax = zeros(1,size_t_plot);
tcounter = 1;
for tcount = floor(linspace(1,prediction_ratio*Nt,size_t_plot))
    Ax(tcounter) = subplot(3,4,tcounter);hold all
    nn=20;
    surface('EdgeColor','red','FaceColor','none',...
        'CData', ones(size(floor(linspace(1,Ny,nn)),2), size(floor(linspace(1,Nx,nn)),2) ),...
        'ZData', ones(size(floor(linspace(1,Ny,nn)),2), size(floor(linspace(1,Nx,nn)),2) ),...
        'YData',y_morphing_fine(floor(linspace(1,Ny,nn)),floor(linspace(1,Nx,nn)),tcount),...
        'XData',x_morphing_fine(floor(linspace(1,Ny,nn)),floor(linspace(1,Nx,nn)),tcount)); axis equal
    nn=60;
%     grey_factor=0.75;
%     surface('EdgeColor',grey_factor*[1,1,1],'FaceColor','none',...'LineStyle','--',...
%         'CData', ones(size(floor(linspace(1,Ny,nn)),2), size(floor(linspace(1,Nx,nn)),2) )-1,...
%         'ZData', ones(size(floor(linspace(1,Ny,nn)),2), size(floor(linspace(1,Nx,nn)),2) )-1,...
%         'YData',y_morphing_fine(floor(linspace(1,Ny,nn)),floor(linspace(1,Nx,nn)),tcount),...
%         'XData',x_morphing_fine(floor(linspace(1,Ny,nn)),floor(linspace(1,Nx,nn)),tcount)); axis equal
    axis off
    title( ['t: ',num2str(tcount),'/',num2str(Nt)] )
    xlabel('x')
    ylabel('y')
    tcounter=tcounter+1;
end
%%
% %%
% figure()
% 
% fields = fieldnames(M_XY);
% for varcount = 1:numel(fields)
%     subplot(2,2,varcount)
%     plot(V_prediction.(fields{varcount})')
% end