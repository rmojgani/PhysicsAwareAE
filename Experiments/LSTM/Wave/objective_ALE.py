#!/usr/bin/env python3
"""
Created on Thu Feb 20 11:21:58 2020

@author: rmojgani
based on: objective_ALE_clean
"""
import numpy as np
from scipy.interpolate import interp1d as interp1
from numpy.linalg import norm

#from scipy.interpolate import interp2d as interp2

def objective_ALE_clean(q,
                        x_0, M_x_0, Nx, M_t_0, Nt, bc,
                        k, X_E_FOM,  
                        Gamma_x , Gamma_t, Gamma_0,
                        Flag, 
                        size_x, ix, ix_fine, 
                        size_t, it, it_fine,
                        boundary_function_map):
    # Decode number of bases
    k_x = k["x"]
    k_w = k["w"]
    # Fine Grid Reconstruction
    [M_x_ALE, dx, U_xx, V_xx] = gridgen_coarse_UV_fine(q, 
                                                       k_x, Flag, Nx, Nt,
                                                       size_x, ix, ix_fine, 
                                                       size_t, it, it_fine,
                                                       M_x_0,
                                                       boundary_function_map);
    # Interpolate the HFM snapshot from Eulerian grid to the ALE grid
    M_w_ALE = interp1d_loop(M_x_0, X_E_FOM, M_x_ALE, Nx, Nt, bc);
    # Project the ALE snapshot 
    M_w_Proj = SVD_project(M_w_ALE, k_w);
    # Interpolate the ALE snapshot from the ALE grid to the Eulerian grid
    X_ROM_Proj = interp1d_loop(M_x_ALE, M_w_Proj, M_x_0, Nx, Nt, bc);
    #
    f_1 = X_E_FOM - X_ROM_Proj; # w on Eul
    f_2 = Gamma_0["x"] * np.dot( Gamma_x, U_xx   );
    f_3 = Gamma_0["t"] * np.dot( Gamma_t, V_xx.T );
    f = np.concatenate(	[f_1.reshape(Nx*Nt), f_2.reshape(Nx*k_x), f_3.reshape(k_x*Nt)] )
    return f#, M_x_ALE, X_E_FOM, X_ROM_Proj, U_d, V_d

def constraint_fun( q, 
                   k_x, Flag, Nx, Nt,
                   size_x, ix, ix_fine, 
                   size_t, it, it_fine,
                   M_x_0,
                   boundary_function_map):
    [M_x_ALE, dx, U_xx, V_xx] = gridgen_coarse_UV_fine(q, 
                               k_x, Flag, Nx, Nt,
                               size_x, ix, ix_fine, 
                               size_t, it, it_fine,
                               M_x_0,
                               boundary_function_map);
    
    diffx = grid_spacing(M_x_ALE)
    return diffx
#%% 
def q2UV( q, Nx, Nt, k, Flag):
    # Extracts and reshape basis functions vectorized in vector q based on its
    # size
    U_d = q[         0:Flag[0]       ].reshape(Nx ,k);
    V_d = q[ Flag[0]:sum(Flag[0:2])  ].reshape(k,Nt );
    
    d = np.dot(U_d,V_d)
    
    return U_d, V_d, d

def USV_reduce(U, S, Vh, k):
    V = np.dot( S[0:k,0:k], Vh[0:k,:] );
    U = U[:  ,0:k]
    return U, V


def SVD_project(M, k, output_UV = False):
    [U, S, Vh] = np.linalg.svd(M, full_matrices=False);
    [U, V] = USV_reduce(U, np.diag(S), Vh, k);
    if output_UV:
        return np.dot( U, V), U, V
    else:
        return np.dot( U, V)

def gridgen_coarse_UV_fine( q, 
                           k_x, Flag, Nx, Nt,
                           size_x, ix, ix_fine, 
                           size_t, it, it_fine,
                           M_x_0,
                           boundary_function_map,
                           mykind='linear'):
    [U_x_coarse , V_x_coarse, d_coarse] = q2UV( q, size_x, size_t, k_x, Flag)
    [U_x_coarse , V_x_coarse] = UV_init_boundary(U_x_coarse , V_x_coarse, boundary_function_map);
    # Construction of fine grid from coarse grid bases
    U_xx = np.zeros((Nx,k_x));
    V_xx = np.zeros((k_x,Nt));
    for kcount in range(k_x):
        U_x_fun = interp1(ix, U_x_coarse[:,kcount], kind = mykind);#cubic, linear, quadratic
        U_xx[:,kcount] = U_x_fun(ix_fine);
        
        V_x_fun = interp1(it, V_x_coarse[kcount,:], kind = mykind);
        V_xx[kcount,:] = V_x_fun(it_fine);
    # Construction of fine grid deviation from Eulerian grid
    dx = np.dot(U_xx, V_xx);
    # Construction of fine grid
    M_x_ALE = dx + M_x_0;
    
    return M_x_ALE, dx, U_xx, V_xx
    
def grid_spacing(M_x):
    return np.diff(M_x.T)

def interp1d_loop(M_x_in, X_in, M_x_Out, Nx, Nt, bc, mykind='linear'):
    X_out = np.zeros((Nx,Nt));
    for tcount in range(Nt):
        state_fun = interp1(M_x_in[:,tcount], X_in[:,tcount], 
                            kind = mykind, fill_value=(bc,bc), bounds_error=False );
        X_out[:,tcount] = state_fun( M_x_Out[:,tcount])
    
    #    state_fun = interp2(M_x_0, M_t_0, X_E_FOM, 
    #                        kind = 'linear', fill_value=(bc,bc), bounds_error=False );
    #    M_w_ALE = state_fun( M_x_ALE, M_t_0)
    return X_out



def UV_init_boundary(U, V, boundary_function_map):
    V = boundary_function_map["UV_init"](V)
    U = boundary_function_map["UV_boundary"](U)
    return U, V
    
def UV_init(V):
     V[:,0] = 0;
     return V
 
def UV_boundary(U):
    U[ 0,:] = 0;
    U[-1,:] = 0;
    return U

#%% Optimization function callbacks
#def callbackF(Xi):
#    global ps
#    ps.append(Xi)
#    print('Norm of current time step:', str(norm(Xi)) )