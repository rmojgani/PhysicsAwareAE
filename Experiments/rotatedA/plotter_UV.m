% SVD Orthogonalization
d_x_pos = bsxfun(@plus, x, -reshape(x_e,[Nx*Ny,1]) );
d_y_pos = bsxfun(@plus, y, -reshape(y_e,[Nx*Ny,1]) );


% d_x_pos = U_x_pos*V_x_pos;
% d_y_pos = U_y_pos*V_y_pos;


[ux,sx,vx]=svd(d_x_pos,'econ');
[uy,sy,vy]=svd(d_y_pos,'econ');

[U_x_pos_ortho, V_x_pos_ortho] = USV_reduce(ux, sx ,vx, n_pos);
[U_y_pos_ortho, V_y_pos_ortho] = USV_reduce(uy, sy ,vy, n_pos);
%% Gram-Schmidt Orthogonalization
% U_x_pos_ortho = gramschmidt_ortho(U_x_pos')';
% U_y_pos_ortho = gramschmidt_ortho(U_y_pos')';
% V_x_pos_ortho = gramschmidt_ortho(V_x_pos')';
% V_y_pos_ortho = gramschmidt_ortho(V_y_pos')';
%% Test Orthogonalization
U_prod = zeros(n_pos, n_pos);
V_prod = zeros(n_pos, n_pos);
for icount=1:1:n_pos
    for jcount=1:1:n_pos
        U_prod(icount, jcount) = U_x_pos_ortho(:,icount)'*U_x_pos_ortho(:,jcount) ;
        V_prod(icount, jcount) = V_x_pos_ortho(icount,:) *V_x_pos_ortho(jcount,:)'; % Note: V_x and V_y are not normalized
    end
end
%% Plot Original Bases
figure();
ix = 2;
iy = 4;
subplot(ix,iy,+1); plot(angleM, V_x_pos);
xlabel(  '$\theta$','Interpreter','latex','FontSize', 19); xlim([angleM(1) angleM(end)])
ylabel(  '$V_{x_{pos}}$','Interpreter','latex','FontSize', 19)

subplot(ix,iy,+2); plot(angleM, V_y_pos);
xlabel(  '$\theta$','Interpreter','latex','FontSize', 19); xlim([angleM(1) angleM(end)])
ylabel(  '$V_{y_{pos}}$','Interpreter','latex','FontSize', 19)

subplot(ix,iy,+3); plot(angleM, V);
xlabel(  '$\theta$','Interpreter','latex','FontSize', 19); xlim([angleM(1) angleM(end)])
ylabel('V_{w} on Eulerian grid');% xlim([1 31])
 
subplot(ix,iy,+4); plot(angleM, V_lag);
xlabel(  '$\theta$','Interpreter','latex','FontSize', 19); xlim([angleM(1) angleM(end)])
ylabel('V_{w} on Morphying grid');%xlim([1 31])

%% Plot Orthogonalized Bases
% figure();
ix = 2;
% iy = 4;
subplot(ix,iy,iy+1); plot(angleM, -V_x_pos_ortho');
xlabel(  '$\theta$','Interpreter','latex','FontSize', 19); xlim([angleM(1) angleM(end)])
ylabel(  '$V_{x_{pos}}$','Interpreter','latex','FontSize', 19)
title(  'Orthogonalized','Interpreter','latex','FontSize', 19);

subplot(ix,iy,iy+2); plot(angleM, -V_y_pos_ortho'); ylabel(  '$V_{y_{pos}}$','Interpreter','latex','FontSize', 19)
xlabel(  '$\theta$','Interpreter','latex','FontSize', 19); xlim([angleM(1) angleM(end)])
title(  'Orthogonalized','Interpreter','latex','FontSize', 19);

% subplot(ix,iy,+3); plot(V'); ylabel('V_{w} on Eulerian grid');% xlim([1 31])
% subplot(ix,iy,+4); plot(V_lag(:,1:n_mag)); ylabel('V_{w} on Morphying grid');%xlim([1 31])