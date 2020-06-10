# -----------------------------------------------------------------------------------------------
import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
import scipy.sparse.linalg as sla
#import numpy.linalg as la
import time as time
import matplotlib.pyplot as plt
from scipy import interpolate
# -----------------------------------------------------------------------------------------------
def weights(z,x,m):
    #  Calculates FD weights. The parameters are:
    #   z   location where approximations are to be accurate,
    #   x   vector with x-coordinates for grid points,
    #   m   highest derivative that we want to find weights for
    #   c   array size m+1,lentgh(x) containing (as output) in
    #       successive rows the0weights for derivatives 0,1,...,m.

    n = len(x)
    c=np.zeros((m+1,n))
    c1=1; c4=x[0]-z; c[0,0]=1;
    for i in range(1,n):
        mn=min(i,m)
        c2=1
        c5=c4
        c4=x[i]-z
        for j in range(0,i):
            c3=x[i]-x[j]
            c2 = c2*c3
            if j == i-1 :
                c[1:mn+1][:,i]= c1*( np.asarray(np.arange(1,mn+1)) * c[0:mn,i-1]-c5*c[1:mn+1][:,i-1] )/c2
                c[0,i] =-c1*c5*c[0,i-1]/c2
            c[1:mn+1][:,j]=(c4*c[1:mn+1][:,j]- np.asarray(np.arange(1,mn+1)) *c[0:mn][:,j] )/c3
            c[0,j]=c4* c[0,j]/c3
        c1=c2;
    return(c)
    
    
def fd_normal(N,order,x,derivative):
    #  Calculates FD derivative matrix. The parameters are:
    #   N           mesh size
    #   order+1     ODD only!!! (3,5,7,...)
    #   x           mesh points
    #   derivative  ^th derivative (1,2,...)
    
    Dm = np.zeros((N, N))
    oo = int(int(order-1)/2)
    for n in np.arange(0,oo, 1):
        w=weights(x[n],x[0:order],derivative)
        Dm[n,0:order]=w[derivative,:]

    for n in np.arange(oo,N-oo, 1):
        w=weights(x[n],x[n-oo:n+oo+1 ],derivative);
        Dm[n,n-oo:n+oo+1]=w[derivative,:]

    for n in np.arange(N-oo,N, 1):
        w=weights(x[n],x[N-order:N+1],derivative)
        Dm[n,N-order:N]=w[derivative,:]
                            
    return (Dm)
    
#%% -----------------------------------------------------------------------------------------------
def snapshot_gen_wave(xmax, tmax, x_mid_0, Nx, Nt, wave_speed, batch_repeat):
    #X_E_FOM = wave_2nd_Eulerian_implicit_2ndorder(Dxx, w_0, Nx, Nt, dt, c, sigma);

    print('Snapshot generator Burgers equation')
    #    xmax = 1.5*3
    #    tmax = 1
    x_snapshot = np.zeros((Nx,Nt))
    y_snapshot = np.zeros((Nx,Nt))
    #    Nx = 500
    #    Nt = 250
    sigma  = 0;
    c = 1.0;
    
    x = np.linspace(0,xmax,Nx)
    dx = x[1]-x[0];
    dt = tmax /(Nt-1.0);    
    print(dx)
    print(dt)
    
    x_grid  = np.arange(0,xmax+dx,dx)
            
    Dxx = fd_normal(Nx,3,x_grid,2);
    
    Dxx  = Dxx[1:-1:1][:,1:-1:1];
    
    Nx     = Nx - 2;
    x_grid = x_grid[1:-1:1];
    
    sI = np.eye(Nx, Nx);

    # To sparse
    Dxx = csc_matrix(Dxx)
    sI  = csc_matrix(sI)
    
    #
    w_0     = np.exp(-((x_grid-x_mid_0)/0.1)**2)
    
    #w       = np.zeros((1,Nx))[0];

    w_back  = w_0.copy();
    w_back2 = w_0.copy();
    w_back3 = w_0.copy();
    
    x_snapshot[1:-1:1,0] = w_0.reshape(Nx,)
    y_snapshot[1:-1:1,0] = w_0.reshape(Nx,)
    
    mytime = 0;

    t = time.time()
    for n in range(1,Nt,1):
        A = (2.0+1.5*sigma*dt) * sI - (c*c) * (dt*dt) * Dxx;
        b = (5.0+2.0*sigma*dt) * w_back - (4.0+0.5*sigma*dt) * w_back2 + w_back3;
        # Sparse Solver 
        w = sla.spsolve(A, b, 'NATURAL') # w = sla.inv(A).dot(b);#
        mytime = mytime+dt;
        #print('n',mytime)
        # Save
        x_snapshot[1:-1:1,n] = w_0.reshape(Nx,)
        y_snapshot[1:-1:1,n] = w.reshape(Nx,)
        # Update
        w_back3 = w_back2.copy();
        w_back2 = w_back.copy();
        w_back  = w.copy();
    
    t = time.time() - t
    print('simulation time:',t)
    print('ttt',mytime)
    
    if wave_speed!=0:
        for n in range(0,Nt,1):
            x_moved = np.linspace(n*dt*wave_speed,xmax+n*dt*wave_speed,Nx+2)
            y = y_snapshot[:,n] 
            f = interpolate.interp1d(x, y, fill_value="extrapolate")
            y_snapshot[:,n]  = f(x_moved)   # use interpolation function returned by `interp1d`

    fig = plt.figure(figsize=(5, 5))
    plt.plot(y_snapshot[:,np.linspace(1,Nt,10).astype(int)-1])
    plt.ylim(-1, 1)
    plt.show()

    x_train = np.zeros((1,Nx+2,Nt))
    y_train = np.zeros((1,Nx+2,Nt))
    #    for batch in range(batch_repeat-):
    #        print('hello')
    batch = 0 
    x_train[batch,:,:] = x_snapshot
    y_train[batch,:,:] = y_snapshot
    
    return x_train, y_train, x, t
#%% -----------------------------------------------------------------------------------------------
def snapshot_decode(xmax, tmax, x_mid_0, Nx, Nt, wave_speed, batch_repeat, y_snapshot, bc):
    print('Snapshot decoder')
    
    print(wave_speed)
    x = np.linspace(0,xmax,Nx)
    t = np.linspace(0,tmax,Nt)
    dx = x[1]-x[0];
    dt = t[1]-t[0];
    
    #    dx = xmax/(Nx-1.0)
    #    dt = tmax/(Nt-1.0)
    
    x_grid  = np.arange(0,xmax+1e-10,dx)
    
#    y_snapshot = np.zeros((Nx+2,Nt))+bc

    if wave_speed!=0:
        for n in range(0,Nt,1):
            x_moved = np.linspace(n*dt*wave_speed,xmax+n*dt*wave_speed,Nx)
            y = y_snapshot[:,n] 
#            f = interpolate.interp1d(x_moved, y, fill_value="extrapolate")
            f = interpolate.interp1d(x_moved, y, kind="previous", fill_value=(bc, bc), bounds_error=False)
            y_snapshot[:,n]  = f(x)   # use interpolation function returned by `interp1d`
        
    return y_snapshot
##%%
### =============================================================================
### =======   #https://www.tensorflow.org/guide/keras/rnn
## tf.keras.layers.SimpleRNN, a fully-connected RNN where the output from 
##                           previous timestep is to be fed to next timestep.
### =============================================================================
#
#epoch_max   = 1001*2*5#*10
#
#Nx = 500+1;
#Ny = 500;
#Nt = 2000+1;
#Encoder_size = 10#*25;
#RNN_size     = 5#*50;
#
#depth = 1
#reg   = 0#1e-9#1e-4
#hyperparams = [RNN_size, depth, reg]
##neurons , depth , reg = int(hyperparams[0]), int(hyperparams[1]), hyperparams[2]
#act_RNN , act_Dense = 'tanh', 'linear'#'tanh', 'linear'
#
##Nx = 50;
##Ny = 50;
##Nt = 25;
##RNN_size = 20;
#
#tmax = 1.0;
#xmax = 1.0;
#x_mid_0 = 0.5;
#batch_repeat = 2;
#wave_speed = 0;
#(x_train_org, y_train_org, x, t) = snapshot_gen_wave(xmax, tmax, x_mid_0, Nx, Nt, wave_speed, batch_repeat)