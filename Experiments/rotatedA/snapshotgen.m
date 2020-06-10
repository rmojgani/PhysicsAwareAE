%% Rotate
if strcmp(my_case,'rotate')
    file_name='A50C';
%     file_name='A25C';
    % file_name='o50';
    
    Img = imread([file_name,'.png']);
    
    nA = length(Img);
    x = 1:1:nA;
    y = 1:1:nA;
    [X,Y] = meshgrid(x,y);
    
    
    Img = im2double(Img(:,:,1))*1.0;
    %     Imgs = imresize(Img, 4); % not any pros
    
    angle.delta = 3 ;
    angle.max   = 90;
    angleM = 0:angle.delta:angle.max;
    
    %     figure();subplot(1,2,1);subimage(Img);subplot(1,2,2);subimage(double(Img))
    
    A = zeros(nA,nA,length(angleM));
    
    if  strcmp(imrotate_method,'griddata')
        A = imrotate_via_griddata(Img, angleM, griddata_method);
    else
        for icount = 1:1:length(angleM)
            %         A(:,:,icount) = 1-fast_rotate(1-Img(:,:,1),angleM(icount));
            
            A(:,:,icount) =  1-imrotate(1-Img, angleM(icount), imrotate_method, 'crop');
            
            %         Imgr = 1-imrotate(1-Imgs,angleM(icount),imrotate_method,'crop');
            %         A(:,:,icount) = imresize(Imgr, 0.25);
        end
    end
elseif strcmp(my_case,'rotateGD')
    file_name='A50C';
    
    Img = imread([file_name,'.png']);
    
    nA = length(Img);
    x = 1:1:nA;
    y = 1:1:nA;
    [X,Y] = meshgrid(x,y);
    
    
    Img = im2double(Img(:,:,1))*1.0;
    
    angle.delta = 3 ;
    angle.max   = 90;
    angleM = 0:angle.delta:angle.max;
    
    [Nx, Ny] = size(Img);    % snapshot size
    Nt = length(angleM);
    x_e = X; %clear X;     % eulerian grid
    y_e = Y; %clear Y;     % eulerian grid
    %%
    size_x = Nx;%5;
    size_y = Ny;%5;
    size_t = Nt;%4;
    %
    angleM = 0:angle.delta:angle.max;
    for icount = 1:1:length(angleM)
        [Xr, Yr] = rotate2Dmesh(X,Y,-angleM(icount));
        xx(:,icount) = Xr(:);
        yy(:,icount) = Yr(:);
    end
    
    d_x = bsxfun(@plus, xx, -reshape(x_e,[Nx*Ny,1]) );
    d_y = bsxfun(@plus, yy, -reshape(y_e,[Nx*Ny,1]) );
    
    x = bsxfun(@plus, d_x, reshape(x_e,[Nx*Ny,1]) );
    y = bsxfun(@plus, d_y, reshape(y_e,[Nx*Ny,1]) );
    %
    M_tilde_j = zeros(Nx*Ny , 1); % Pre-allocation
    M_tilde_lag = zeros(Nx*Ny , Nt); % Pre-allocation
    for j = 1:Nt
        xx = reshape(x(:,j),[Nx,Ny]);
        yy = reshape(y(:,j),[Nx,Ny]);
        M_tilde_j = griddata_fill(x_e, y_e, Img, xx, yy, griddata_method);
        M_tilde_lag(:,j) = M_tilde_j(:);
    end
    % Interpolation
    M_tilde_lag = bsxfun(@times, M_tilde_lag(:,1), ones(Nx*Ny,Nt) );
    %
    M_tilde_lag_back = zeros(Nx, Ny , Nt); % Pre-allocation
    for j = 1:Nt
        xx = reshape(x(:,j),[Nx,Ny]);
        yy = reshape(y(:,j),[Nx,Ny]);
        ww_lag = reshape(M_tilde_lag(:,j),[Nx,Ny]);
        M_tilde_lag_back(:,:,j)  = griddata_fill(xx, yy, ww_lag, ...
            x_e, y_e, griddata_method);
    end
    
    M = reshape(M_tilde_lag_back, [Nx,Ny, Nt]);
elseif strcmp(my_case,'RotMorph')
    %% Rotate and Morph
    %     file_name = [my_case,'0O250'];
    %     Img1 = imread('0250.png');
    %     Img2 = imread('O250.png');
    file_name = [my_case,'0A50'];
    Img1 = imread('050.png');
    Img2 = imread('A50C.png');
    
    nA = length(Img1);
    x = 1:1:nA;
    y = 1:1:nA;
    [X,Y] = meshgrid(x,y);
    
    Img1 = im2double(Img1(:,:,1))*1;
    Img2 = im2double(Img2(:,:,1))*1;
    
    angle.delta = 15*4;
    angle.max   = 90*4;
    angleM = 0:angle.delta:angle.max;
    
    A = zeros(nA,nA,length(angleM));
    for icount = 1:1:length(angleM)
        A_morph = (icount-1)/(length(angleM)-1) * Img1 +  (length(angleM) - icount +1)/(length(angleM)-1)*Img2;
        A(:,:,icount) = 1-imrotate(1-A_morph,angleM(icount),imrotate_method,'crop');
    end
end
M=A;