% clear
% load('main_solve_sm_case12.mat')
% load('main_solve_sm_k=8_Proj.mat')
% load('main_solve_sm_k=2_ROM.mat')
% load('main_solve_sm_k=12_ROM')
% kM = 14;%:2:10%:2:12;
% load('main_solve.mat')
%%
fig1 = figure();
fields = fieldnames(error.Eul_Proj);
for filecount = 1:1:numel(fields)
    subplot(1,numel(fields),filecount)
    title(fields{filecount})
    hold on
    loglog(kM,  error.Eul_Proj.(fields{filecount}), '.-k', 'DisplayName','$\| U_{w_{Eul}} V_{w_{Eul}} - {w}_{Eul} \|$')
    loglog(kM,  error.Eul_ROM.(fields{filecount}), 'o--k', 'DisplayName','$\| U_{w_{Eul}} a_{w_{Eul}} - {w}_{Eul} \|$')
    loglog(kM,  error.ALE_Proj_PALE_2_Eul.(fields{filecount}), 'o-m', 'DisplayName','$\| P(U_{w_{ALE}}) V_{w_{ALE}} - {w}_{Eul} \|$')
    loglog(kM,  error.ALE_ROM_PALE_2_Eul.(fields{filecount}), 'o--c', 'DisplayName','$\| P(U_{w_{ALE}}) a_{w_{ALE}} - {w}_{Eul} \|$')
    
%     loglog(kM, error.ALE_Proj_ALE2ALE, '.-r', 'DisplayName','$\| U_{w_{ALE}} V_{w_{ALE}} - {w}_{ALE} \|$')
%     loglog(kM, error.ALE_Proj_PALE_2_PALE, '*-r', 'DisplayName','$\| P(U_{w_{ALE}}) V_{w_{ALE}} - P({w}_{ALE}) \|$')
%     loglog(kM, error.ALE_ROM_PALE_2_PALE, '*--g', 'DisplayName','$\| P(U_{w_{ALE}}) a_{w_{ALE}} - P({w}_{ALE}) \|$')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    lgd = legend('Location','northeast');
    set(lgd,'Interpreter','latex','FontSize',14);
    xlabel('$k_w$','Interpreter','latex','FontSize',20);
    ylabel('$error$','Interpreter','latex','FontSize',20);
    % xlim([0 10])
    grid on
end
pause(5.005)
%%
% fields = fieldnames(M_Eulerian);
nnt = size(M_Eulerian.p,2);
for filecount = 1:numel(fields)
    M_Eulerian_XY.(fields{filecount}) = reshape(  M_Eulerian.(fields{filecount}), [Nx,Ny,nnt]);
    M_Eulerian_Proj_XY.(fields{filecount}) = reshape(  M_Eulerian_Proj.(fields{filecount}), [Nx,Ny,nnt]);
    M_ALE_Proj_on_stationary_grid_XY.(fields{filecount}) = reshape(  M_ALE_Proj_on_stationary_grid.(fields{filecount}), [Nx,Ny,nnt]);
end
for filecount = 1:numel(fields)
    M_Eulerian_XY.(fields{filecount}) = reshape(  M_Eulerian.(fields{filecount}), [Nx,Ny,nnt]);
    M_Eulerian_ROM_XY.(fields{filecount}) = reshape(  M_Eulerian_ROM.(fields{filecount}), [Nx,Ny,nnt]);
    M_ALE_ROM_on_stationary_grid_XY.(fields{filecount}) = reshape(  M_ALE_ROM_on_stationary_grid.(fields{filecount}), [Nx,Ny,nnt]);
end
%%
size_t_plot = 10;
fig2 = figure();
Ax = zeros(1,size_t_plot);
tcouner = 1;
filecount = 3;
for tcount = floor(linspace(1,Nt,size_t_plot))
    Ax(tcouner) = subplot(5,2,tcouner);
    plot(diag(M_Eulerian_XY.(fields{filecount})(:,:,tcount))                    ,' -k'); hold on
    plot(diag(M_ALE_Proj_on_stationary_grid_XY.(fields{filecount})(:,:,tcount)) ,'.:r');
    plot(diag(M_Eulerian_Proj_XY.(fields{filecount})(:,:,tcount))               ,'.-b');
%     hold off
    drawnow
%     pause(0.001)
    tcouner=tcouner+1;
    ylabel(fields{filecount})
end
xlabel('Nx')
linkaxes(Ax,'xy')
% ylim([-1 1])
%% M_Eulerian_ROM_XY
tcount = Nt;
fig3 = figure();
for varcount = 1:numel(fields)
    subplot(3, numel(fields), varcount);
    surf(X, Y, M_Eulerian_XY.(fields{varcount})(:,:,tcount) )
    axis square, shading flat
    view([0 90])
    title( [(fields{varcount}),', t:',num2str(tcount),'/',num2str(Nt)] )
    
    subplot(3, numel(fields), numel(fields)+varcount);
    surf(X, Y, M_Eulerian_Proj_XY.(fields{varcount})(:,:,tcount) )
    axis square, shading flat
    view([0 90])
    title( [(fields{varcount}),', t:',num2str(tcount),'/',num2str(Nt)] )

    subplot(3, numel(fields), 2*numel(fields)+varcount);
    surf(X, Y, M_ALE_Proj_on_stationary_grid_XY.(fields{varcount})(:,:,tcount) )
    axis square, shading flat
    view([0 90])
    title( [(fields{varcount}),', t:',num2str(tcount),'/',num2str(Nt)] )
    
end
%%
saveas(fig1,[num2str(1),'.png'])
saveas(fig2,[num2str(2),'.png'])
saveas(fig3,[num2str(3),'.png'])
%%
tcount = nnt;
fig3 = figure();
colorbar_show='on';% 'on' / 'off'
for varcount = 1:numel(fields)
    top1 = -10;
    bottom1 = 10;

    subplot(3, numel(fields), varcount);
    surf(X, Y, M_Eulerian_XY.(fields{varcount})(:,:,tcount) )
    axis square, shading flat
    view([0 90])
    title( [(fields{varcount}),', t: ',num2str(tcount),'/',num2str(Nt)] )
    
    subplot(3, numel(fields), numel(fields)+varcount);
    surf(X, Y, M_Eulerian_ROM_XY.(fields{varcount})(:,:,tcount) )
    axis square, shading flat
    view([0 90])
    title( [(fields{varcount}),', t: ',num2str(tcount),'/',num2str(Nt)] )

    subplot(3, numel(fields), 2*numel(fields)+varcount);
    surf(X, Y, M_ALE_ROM_on_stationary_grid_XY.(fields{varcount})(:,:,tcount) )
    axis square, shading flat
    view([0 90])
    title( [(fields{varcount}),', t: ',num2str(tcount),'/',num2str(Nt)] )
    
    
    bottom    = min( min( min(min(min(M_Eulerian_XY.(fields{varcount})))),...
                min(min(min(M_Eulerian_XY.(fields{varcount})))) ),...
                min(min(min(M_Eulerian_XY.(fields{varcount})))) );
    bottom = min(bottom1, bottom);
    
    top    = max( max( max(max(max(M_Eulerian_XY.(fields{varcount})))),...
                max(max(max(M_Eulerian_XY.(fields{varcount})))) ),...
                max(max(max(M_Eulerian_XY.(fields{varcount})))) );
    top = max(top1, top);

    subplot(3, numel(fields), varcount);
    colorbar('Visible',colorbar_show,'Limits',[bottom top]);

    subplot(3, numel(fields), numel(fields)+varcount);
    colorbar('Visible',colorbar_show,'Limits',[bottom top]);
    
    subplot(3, numel(fields), 2*numel(fields)+varcount);
    colorbar('Visible',colorbar_show,'Limits',[bottom top]);
end
%%
size_t_plot = 10;
fig2 = figure();
Ax = zeros(1,size_t_plot);
tcouner = 1;
varcount = 1;

top1 = -10;
bottom1 = 10;
tcount = Nt;

% bottom   = min( min( min(min(min(M_Eulerian_XY.(fields{varcount})))),...
%             min(min(min(M_Eulerian_ROM_XY.(fields{varcount})))) ),...
%             min(min(min(M_Eulerian_ROM_XY.(fields{varcount})))) );
% bottom = min(bottom1, bottom);
% top    = max( max( max(max(max(M_Eulerian_XY.(fields{varcount})))),...
%             max(max(max(M_Eulerian_ROM_XY.(fields{varcount})))) ),...
%             max(max(max(M_Eulerian_ROM_XY.(fields{varcount})))) );
% top = max(top1, top);

for tcount = floor(linspace(1,Nt,size_t_plot))
    Ax(tcouner) = subplot(5,2,tcouner);
    plot(diag(M_Eulerian_XY.(fields{varcount})(:,:,tcount))                   ,' -k'); hold on
    plot(diag(M_ALE_ROM_on_stationary_grid_XY.(fields{varcount})(:,:,tcount)) ,'.:r');
    plot(diag(M_Eulerian_ROM_XY.(fields{varcount})(:,:,tcount))               ,'.-b');
%     hold off
    drawnow
    pause(0.1)
    tcouner=tcouner+1;
    ylabel(fields{varcount})
%     ylim([floor(bottom*10)/10 ceil(top*10)/10])
    xlim([0 Nx])

end
xlabel('Nx')
% linkaxes(Ax,'xy')