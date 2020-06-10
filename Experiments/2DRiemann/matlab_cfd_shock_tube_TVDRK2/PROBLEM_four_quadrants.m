gamma = 1.4;

% Nx  = 100;
% Ny  = 100;

% CFL = 0.4;

% write_snapshots = false;
% WRITE_DIRECTORY = '/home/maciej/data/';

x = 0:xmax/(Nx-1):xmax;
y = 0:ymax/(Ny-1):ymax;

[X,Y] = meshgrid(x_fine,y_fine);
X=X';
Y=Y';

dx = diff(X,1,1); dx(Nx,:) = dx(Nx-1,:);
dy = diff(Y,1,2); dy(:,Ny) = dy(:,Ny-1);


            % |------------ left ---------|    |----------- right -----------| 
            % p     rho     u       v           p       rho     u       v
dd{1} =     [1      1       0       0           1/10    1/8     0       0;
            1       1       0       0           1/10    1/8     0       0]; 

            % |------------ left ---------|    |----------- right -----------| 
            % p     rho     u       v           p       rho     u       v
dd{3} =     [0.3    0.5323  1.206       0       1.5     1.5     0       0;
            0.029   0.138   1.206       1.206   0.3     0.5323  0       1.206];
        
            % |------------ left ---------|    |----------- right -----------| 
            % p     rho     u       v           p       rho     u       v
dd{6} =     [1.0    2.0     0.75    0.5         1.0     1.0     0.75   -0.5 ;
            1.0     1.0    -0.75    0.5         1.0     3.0    -0.75   -0.5]; 
        
            % |------------ left ---------|    |----------- right -----------| 
            % p     rho     u       v           p       rho     u       v
dd{11} =    [0.4    0.5313  0.8276  0           1.0     1.0     0.1     0;
            0.4     0.8     0.1     0           0.4     0.5313  0.1     0.7276];
        
            % |------------ left ---------|    |----------- right -----------| 
            % p     rho     u       v           p       rho     u       v
dd{12} =    [1      1       0.7276  0           0.4     0.5313  0       0;
            1       0.8     0       0           1.0     1.0     0       0.7276];


            % |------------ left ---------|    |----------- right -----------| 
            % p     rho     u       v           p       rho     u       v
dd{15} =    [ 0.4   0.5197  -0.6259 -0.3        1.0     1.0     0.1     -0.3;
            0.4     0.8     0.1     -0.3        0.4     0.5313  0.1     0.4276];
        
        
            % Case 50: https://computation.llnl.gov/projects/blast/triple-point-shock-interaction
            % |------------ left ---------|    |----------- right -----------| 
            % p     rho     u       v           p       rho     u       v
dd{50} =    [ 1     1       0       0           0.1     0.125   0     0 ;
              1     1       0       0           0.1     1.000   0     0];
        
ts{1} = 0.2; ts{3}=0.8; ts{4}=0.25; ts{6}=0.3; ts{11}=0.75; ts{12}=0.30; ts{15}=10*0.2; ts{17}=0.3;
ts{50} = 5.0;

t_max = ts{my_case};

    clear val
    val.x = 0.75;
    val.y = val.x;
    if any(my_case == [6, 12 , 15])
        val.x = 0.5;
        val.y = val.x;
    elseif my_case == 50;
        val.x = 1/6*xmax;
        val.y = 0.5*ymax;
    end
    m(:,:,1) = and(X< val.x,Y>=val.y); %upper left
    m(:,:,2) = and(X>=val.x,Y>=val.y); %upper right
    m(:,:,3) = and(X< val.x,Y< val.y); %lower left
    m(:,:,4) = and(X>=val.x,Y< val.y); %lower rigth

    p = zeros(Nx,Ny);
  rho = zeros(Nx,Ny);
    u = zeros(Nx,Ny);
    v = zeros(Nx,Ny);
    
    % UPPER LEFT
    p(m(:,:,1)) = dd{my_case}(1,1);
  rho(m(:,:,1)) = dd{my_case}(1,2);
    u(m(:,:,1)) = dd{my_case}(1,3);
    v(m(:,:,1)) = dd{my_case}(1,4);

    p(m(:,:,2)) = dd{my_case}(1,5);
  rho(m(:,:,2)) = dd{my_case}(1,6);
    u(m(:,:,2)) = dd{my_case}(1,7);
    v(m(:,:,2)) = dd{my_case}(1,8);

    p(m(:,:,3)) = dd{my_case}(2,1);
  rho(m(:,:,3)) = dd{my_case}(2,2);
    u(m(:,:,3)) = dd{my_case}(2,3);
    v(m(:,:,3)) = dd{my_case}(2,4);

    p(m(:,:,4)) = dd{my_case}(2,5);
  rho(m(:,:,4)) = dd{my_case}(2,6);
    u(m(:,:,4)) = dd{my_case}(2,7);
    v(m(:,:,4)) = dd{my_case}(2,8);

% make conservative variables
N = Nx*Ny;
q(:,:,1) = rho;
q(:,:,2) = rho.*u;
q(:,:,3) = rho.*v;
q(:,:,4) = p./(gamma-1) + 0.5*rho.*(u.^2 + v.^2);

% clean up
clear rho u v p;