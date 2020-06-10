% vim:set tw=0 ft=matlab ts=2 sw=2 sts=2 et :

% clear all
% clc
% addpath matlab_functions
addpath_host()
display('ssssssssssssssssssssssssssssssssssssssssssssss')
global gamma
global dx dy dz
% global Nx Ny Nz
global Dx Dy Dxx Dyy Dxxxx Dyyyy

% default plot attributes
set(0,'defaultfigurecolor',     [1 1 1]);
set(0,'defaultfigureposition',  [10 700 1700 350])

% ======================================================
%
%     x-direction stored in columns
%
%   ---------------------------------> y
%   |
%   | [ M(1,1)  M(1,2) ...  M(1,Ny) ]
%   | [ M(2,1)  M(2,2) ...  M(2,Ny) ]
%   | [                             ]
%   | [   .     .       .      .    ]
%   | [   .     .        .     .    ]
%   | [   .     .         .    .    ]
%   | [                             ]
%   | [ M(Nx,1) M(N,N) ...  M(Nx,Ny)]
%   V
%
%   x
%
% ======================================================

% problems available:
% *_shock_bubble
% *_four_quadrants
% *_triple_point
% *_explosion
% *_sod

PROBLEM_four_quadrants

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
order_D1 = 11;
Dx = sparse(fd_normal(Nx,order_D1,x,1));
D_one = central_2nd_order(x,0);
Dx(1:order_D1,:) = D_one(1:order_D1,:);
Dx(end-order_D1:end,:) = D_one(end-order_D1:end,:);

Dy = sparse(fd_normal(Ny,order_D1,y,1));
D_one = central_2nd_order(y,0);
Dy(1:order_D1,:) = D_one(1:order_D1,:);
Dy(end-order_D1:end,:) = D_one(end-order_D1:end,:);

order_D2 = 3;
Dxx = sparse(fd_normal(Nx,order_D2,x,2));
D_one = central_2nd_order_second_derivative(x,0);
Dxx(1:order_D2,:) = D_one(1:order_D2,:);
Dxx(end-order_D2:end,:) = D_one(end-order_D2:end,:);

Dyy = sparse(fd_normal(Ny,order_D2,y,2));
D_one = central_2nd_order_second_derivative(y,0);
Dyy(1:order_D2,:) = D_one(1:order_D2,:);
Dyy(end-order_D2:end,:) = D_one(end-order_D2:end,:);

order_D4 = 5;
Dxxxx = sparse(fd_normal(Nx,order_D4,x,4));
Dyyyy = sparse(fd_normal(Ny,order_D4,y,4));

%break;
TIME_INT = 4;

figure;
time_in_loop =0; tic;
t(1) = 0; k = 0;


M_XY.rho = zeros(Nx,Ny,t_save_size);
M_XY.u = zeros(Nx,Ny,t_save_size);
M_XY.v = zeros(Nx,Ny,t_save_size);
M_XY.p = zeros(Nx,Ny,t_save_size);

M_XY.rho(:,:,1) = q(:,:,1);
M_XY.u(:,:,1) = q(:,:,2)./q(:,:,1);
M_XY.v(:,:,1) = q(:,:,3)./q(:,:,1);
M_XY.p(:,:,1) = (gamma - 1 ) * ( q(:,:,4) - 0.5 * M_XY.rho(:,:,1) .* ( M_XY.u(:,:,1).^2 + M_XY.v(:,:,1).^2 ) );
t_save=1;
for tcount=2:1:Nt
% tcount=1;
% while t(end) < tmax;
%     cfl_condition; tcount=tcount+1;
    t(end+1) = t(end) + dt;
    
    if t(end) > t_max
        dt = t_max - t(end-1);
        t(end) = t_max;
    end
    
    fprintf('\rt = %.6f/%.2f (%.2f/%.2f min)',t(end),t_max,time_in_loop,(t_max/t(end))*time_in_loop);
    if  TIME_INT == 1
        q = q + dt*euler_ode(q, Nx, Ny);
    elseif TIME_INT == 2
        q1 = q + dt*euler_ode(q, Nx, Ny);
        q  = 0.5*q + 0.5*q1 + 0.5*dt*euler_ode(q1, Nx, Ny);
    elseif TIME_INT == 3
        q1 = q + dt*euler_ode(q, Nx, Ny);
        q2 = (3/4)*q + (1/4)*q1 + (1/4)*dt*euler_ode(q1, Nx, Ny);
        q  = (1/3)*q + (2/3)*q2 + (2/3)*dt*euler_ode(q2, Nx, Ny);
    elseif TIME_INT == 4
        k1 = euler_ode(q, Nx, Ny);
        k2 = euler_ode(q + dt/2*k1, Nx, Ny);
        k3 = euler_ode(q + dt/2*k2, Nx, Ny);
        k4 = euler_ode(q + dt*k3, Nx, Ny);
        q = q + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    end
    
    % Boundary 
    [p, u, v, rho] = Q_2_primitive_2D(q, gamma);
    
    % zero derivative (gradient) bc for p,T,rho
    % zero velocity bc (U(1)+U(2))/2 = 0
    rho( :, 1) = rho(:,   2);
    rho( :,Ny) = rho(:,Ny-1);
    rho( 1, :) = rho(2  , :);
    rho(Nx, :) = rho(Nx-1,:);
 
%     u(: , 1) = -u(:,   2);
%     u(: ,Nx) = -u(:,Nx-1);
%     u( 1, 2:Ny-1) = -u(2   , 2:Ny-1);
%     u(Nx, 2:Ny-1) = -u(Nx-1,2:Ny-1);
%     
%     v(: , 1) = -v(:,   2);
%     v(: ,Nx) = -v(:,Nx-1);
%     v( 1, 2:Ny-1) = -v(2   , 2:Ny-1);
%     v(Nx, 2:Ny-1) = -v(Nx-1,2:Ny-1);
    
    p( :, 1) = p(:,   2);
    p( :,Ny) = p(:,Ny-1);
    p( 1, :) = p(2  , :);
    p(Nx, :) = p(Nx-1,:);
    % Primitve to Flux
    q = primitive_2_Q_2D(p, u, v, rho, gamma);
    
    if or( mod(tcount-1, t_save_period)==0, tcount == Nt )
        t_save = t_save + 1;
        [ M_XY.p(:,:,t_save), M_XY.u(:,:,t_save), M_XY.v(:,:,t_save), M_XY.rho(:,:,t_save)] = Q_2_primitive_2D(q, gamma);
    end
    
%     if k == 5e1;
%         subplot(1,2,1);
% 
%         plot3( [xmax/2 xmax/2], [     0     ymax  ], [1000 1000], 'r')
%         hold on
%         plot3( [   0   xmax  ], [  ymax/3   ymax/3], [1000 1000], '.-g', 'LineWidth',2)
%         plot3( [   0   xmax  ], [  ymax/2   ymax/2], [1000 1000], 'g')
%         plot3( [   0   xmax  ], [2*ymax/3 2*ymax/3], [1000 1000], ':g', 'LineWidth',2)
%         surf(X,Y,q(:,:,1));
%         shading flat; axis equal; colorbar; axis tight;
%         view([0 90])
%         title(['t=',num2str((tcount-1)*dt)])
%         hold off
%         
%         subplot(2,6,4)
%         plot(diag(q(:,:,1)));
%         ylabel('$\rho$','Interpreter','latex')
%         
%         subplot(2,6,6)
%         plot( q(end/2,:,1), (1:1:Ny)*ymax/Ny, 'r')
%         xlabel('$\rho$','Interpreter','latex')
%         ylabel('$y$','Interpreter','latex')
% 
%         subplot(2,6,10)
%         plot((1:1:Nx)*xmax/Nx, q(:,end/3,1) , '.-g', 'LineWidth',2)
%         xlabel('$x$','Interpreter','latex')
%         ylabel('$\rho$','Interpreter','latex')
% 
%         subplot(2,6,11)
%         plot((1:1:Nx)*xmax/Nx, q(:,end/2,1) , 'g')
%         xlabel('$x$','Interpreter','latex')
%         ylabel('$\rho$','Interpreter','latex')
% 
%         subplot(2,6,12)
%         plot((1:1:Nx)*xmax/Nx, q(:,2*end/3,1) , ':g', 'LineWidth',2)
%         xlabel('$x$','Interpreter','latex')
%         ylabel('$\rho$','Interpreter','latex')
% 
%         drawnow;
%         k = 0;
%     end
%     k = k+1;
    time_in_loop = toc/60;
end
toc
%%
clear k1 k2 k3 k4
%%
fig = figure();

subplot(2,2,1);
    surf(X,Y,M_XY.rho(:,:,1));
subplot(2,2,2);
    surf(X,Y,M_XY.rho(:,:,floor(end/4)));
subplot(2,2,3);
    surf(X,Y,M_XY.rho(:,:,floor(end/2)));
subplot(2,2,4);
    surf(X,Y,M_XY.rho(:,:,end));

for i =1:1:4
    subplot(2,2,i)
    shading flat; axis equal; colorbar; axis tight;    view([0 90])
end
saveas(fig,[num2str(Nx),'.png'])
%%
fig = figure();

subplot(2,2,1);
    plot(diag(M_XY.rho(:,:  , floor(end/2) )));
subplot(2,2,2);
    plot(     M_XY.rho(floor(Nx*2/3),:, floor(end/2) ));
subplot(2,2,3);
    plot(diag(M_XY.rho(:,:     ,end)));
subplot(2,2,4);
    plot(     M_XY.rho(floor(Nx*2/3),:,end));
    
saveas(fig,['line_',num2str(Nx),'.png'])