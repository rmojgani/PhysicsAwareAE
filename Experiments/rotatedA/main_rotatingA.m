close all;
clear;
clc;
%% ******************************************************************************
% solving:
%
%    min      |M-P{U_w*V_w}| + |Gamma_x*(U_x*V_x)| + |(U_x*V_x)*Gamma_t'|
%  U_x,V_x    st:   (1) rank(U_w*V_w) = n_mag
%                   (2) rank(U_x*V_x) = n_pos
%
%******************************************************************************
griddata_method = 'linear'; % 'linear', 'nearest', 'natural', 'cubic', or 'v4'
imrotate_method = 'griddata'; % 'nearest', 'bilinear', 'bicubic', or 'griddata'
my_case = 'rotate'; % 'rotate', 'RotMorph'
%% ******************************************************************************
snapshotgen
%% ******************************************************************************
[Nx, Ny, Nt] = size(M);    % snapshot size
x_e = X; clear X;     % eulerian grid
y_e = Y; clear Y;     % eulerian grid

rng(4); %not ok [1-6]

size_x = 7;
size_y = 7;

size_t = Nt;

nA = length(A);
x = 1:1:nA;
y = 1:1:nA;
[X,Y] = meshgrid(x,y);

X = X(floor(linspace(1,Nx,size_x)),floor(linspace(1,Ny,size_y)));
Y = Y(floor(linspace(1,Nx,size_x)),floor(linspace(1,Ny,size_y)));
x_e_c = X;
y_e_c = Y;

n_pos = 2; % number of modes
n_mag = 1; % number of modes
k = [n_pos, n_mag];
%% ******************************************************************************
Gamma_0_x = 100;
Gamma_0_t = 100;
Gamma_x = Gamma_0_x*sparse(fd_normal(size_x,3,x_e_c(1,:),2));
Gamma_t = Gamma_0_t*sparse(fd_normal(Nt,3,[1:Nt],2));

Gamma_x = kron(speye(size_x),Gamma_x) + kron(Gamma_x,speye(size_x));

for icount = 1:1:Nt
    MC(:,icount) = reshape(M(:,:,icount),[Nx*Ny,1]);
end
%%
rng(4);
U_x_coarse = 1*(1-2*rand(size_x*size_y, n_pos ));
V_x_coarse = 1*(1-2*rand(n_pos, Nt    ));

U_y_coarse = 1*(1-2*rand(size_x*size_y, n_pos ));
V_y_coarse = 1*(1-2*rand(n_pos, Nt    ));
%%
q.x    = [ U_x_coarse(:)         ; V_x_coarse(:)         ];
Flag.x = [ length(U_x_coarse(:)) , length(V_x_coarse(:)) ];

q.y    = [ U_y_coarse(:)         ; V_y_coarse(:)         ];
Flag.y = [ length(U_y_coarse(:)) , length(V_y_coarse(:)) ];

q = [q.x; q.y];
%%
delete(gcp('nocreate'))
options = optimoptions(...
    'lsqnonlin',...
    'Algorithm',    'levenberg-marquardt',...
    'Display',      'iter',...
    'MaxFunEvals',  1e6,...
    'MaxIter',      0,...
    'ScaleProblem', 'none',...
    'FinDiffType',  'forward',...
    'UseParallel',  false,...%true,...
    'PlotFcn',      @optimplotresnorm,... % plotting let us stop/pause lsqonlin at any point
    'StepTolerance',1e-12,...
    'TolFun',       1e-12,...
    'InitDamping',  1,...
    'ScaleProblem','jacobian');
%%
F = @(q) objective_rotA_proj2(q,...
    k, Flag, ...%x_coarse, y_coarse, ...
    M, ...
    x_e     , y_e     , ...
    x_e_c   , y_e_c   , ...
    Nx    , Ny    , Nt    ,...
    size_x, size_y, size_t,...
    Gamma_x, Gamma_t, ...
    griddata_method);
options.MaxIter = 1000;              % Maximum number of iterations
tic
q = lsqnonlin(F,q,[],[],options); % added to drop the first iterations from the Residual plot
toc
q = lsqnonlin(F,q,[],[],options);
%% SAVE 
save( [file_name,'MAX',num2str(angle.max),'|',num2str(angle.delta),...
    '_n_mag_',num2str(n_mag),'_n_pos_',num2str(n_pos),'_cp',num2str(size_x),'x',num2str(size_y),...
    ['|G',num2str(Gamma_0_x),'_',num2str(Gamma_0_t)],'|00XX'] )
delete(gcp)