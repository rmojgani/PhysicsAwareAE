# -----------------------------------------------------------------------------------------------
import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix
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
    
def Residual_E(w, w_back, w_0, dt, nu0, Dx, Dx_no_bc, Dxx , Dxx_no_bc):
    R = (w - w_back)/dt + 0.5*Dx.dot(w**2) + Dx.dot(w_0*w) + 0.5*Dx_no_bc.dot(w_0**2) - nu0*Dxx.dot(w) - nu0*Dxx_no_bc.dot(w_0); 
    return(R)
#%% -----------------------------------------------------------------------------------------------
def snapshot_gen_burgers(xmax, tmax, x_mid_0, Nx, Nt, wave_speed, batch_repeat):
    print('Snapshot generator Burgers equation')
    #    xmax = 1.5*3
    #    tmax = 1
    
    #    Nx = 500
    #    Nt = 250
    
    x = np.linspace(0,xmax,Nx)
    t = np.linspace(0,tmax,Nt)
    dx = x[1]-x[0];
    dt = t[1]-t[0];
    
    #    dx = xmax/(Nx-1.0)
    #    dt = tmax/(Nt-1.0)
    
    x_grid  = np.arange(0,xmax+1e-10,dx)
    
    bc = 0.8;
    
    nu0 = 1e-3;
    
    Dx  = fd_normal(Nx,3,x_grid,1);
    Dxx = fd_normal(Nx,3,x_grid,2);
    
    Dx  = Dx[1:-1:1][:,1:-1:1];
    Dxx  = Dxx[1:-1:1][:,1:-1:1];
    
    Nx     = Nx - 2;
    x_grid = x_grid[1:-1:1];
    
    Dx_no_bc  = fd_normal(Nx,3,x_grid,1);
    Dxx_no_bc = fd_normal(Nx,3,x_grid,2);
    
    # To sparse
    Dx=csr_matrix(Dx)
    Dxx=csr_matrix(Dxx)
    Dx_no_bc=csr_matrix(Dx_no_bc)
    Dxx_no_bc=csr_matrix(Dxx_no_bc)
    
    #
    w_0 = 0.5*np.exp(-((x_grid-x_mid_0)/0.1)**2) + bc
    x_0 = x_grid;
    
    
    w      = np.zeros((1,Nx))[0]
    w_back = np.zeros((1,Nx))[0]
    
    x_snapshot = np.zeros((Nx+2,Nt))+bc
    y_snapshot = np.zeros((Nx+2,Nt))+bc
    y_snapshot[1:-1,0] = w_0.copy()
    y_snapshot[1:-1,0] = w_0.copy()
    # Sparse Solver 
    t =time.time()
    for n in range(0,Nt,1):
        for i in range(0,4,1):
            R = Residual_E(w, w_back, w_0, dt, nu0, Dx, Dx_no_bc, Dxx , Dxx_no_bc)
            J = sp.eye(Nx,Nx)/dt + Dx.dot(sp.spdiags(w,0,Nx,Nx)) + Dx.dot(sp.spdiags(w_0,0,Nx,Nx)) - nu0*Dxx
            dw = sla.spsolve(J, R) # w = w - numpy.linalg.inv(J).dot(R)
            w = w - dw  
        w_back = w.copy()
        x_snapshot[1:-1,n] = w_0.copy()
        y_snapshot[1:-1,n] = w_0+w_back.copy()
    t = time.time() - t
    print(t)
    w_back_spsolver = w_back.copy()
    
    if wave_speed!=0:
        for n in range(0,Nt,1):
            x_moved = np.linspace(n*dt*wave_speed,xmax+n*dt*wave_speed,Nx+2)
            y = y_snapshot[:,n] 
            f = interpolate.interp1d(x, y, fill_value="extrapolate")
            y_snapshot[:,n]  = f(x_moved)   # use interpolation function returned by `interp1d`

    fig = plt.figure(figsize=(5, 5))
    plt.plot(y_snapshot[:,0:-1:10])
    plt.ylim(0.7, 1.3)
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
