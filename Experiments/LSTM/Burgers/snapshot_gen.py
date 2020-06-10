#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 13:44:48 2020

@author: rmojgani
"""

import numpy as np

def snapshot_gen_v1(xmax, tmax, x_mid_0, Nx, Nt, wave_speed, batch_repeat):
    x = np.linspace(0,xmax,Nx)
    t = np.linspace(0,tmax,Nt)
    dt = t[2]-t[1];
    x_snapshot = np.zeros((Nx,Nt))
    for tcount in range(Nt):
        x_mid = x_mid_0 + wave_speed *dt * (tcount)        
        x_snapshot[:,tcount] =  (1-0*0.5*dt * (tcount-1)) * np.exp(-2*np.power(x-x_mid,2)/(0.1+0.15)**2)     
    
    x_train = np.zeros((batch_repeat,Nx,Nt))
    for batch in range(batch_repeat-1):
        x_train[batch,:,:] = x_snapshot
    
    y_train = np.copy(x_train)

    return x_train, y_train, x, t


def snapshot_gen(xmax, tmax, x_mid_0, Nx, Nt, wave_speed, batch_repeat):
    x = np.linspace(0,xmax,Nx)
    t = np.linspace(0,tmax,Nt)
    dt = t[1]-t[0];
            
    x_snapshot = np.zeros((Nx,Nt))
    y_snapshot = np.zeros((Nx,Nt))
    for tcount in range(Nt):
        x_mid = x_mid_0 + wave_speed *dt * (tcount)
        #x_mid = x_mid_0 + 0.25*np.sin(10*wave_speed *dt * (tcount))
        x_snapshot[:,tcount] =  (1-0*0.5*dt * (tcount-1)) * np.exp(-2*np.power(x-x_mid_0,2)/0.1**2)     
        y_snapshot[:,tcount] =  (1-0.5*dt * (tcount-1)) * np.exp(-2*np.power(x-x_mid,2)/0.1**2)     
    
    x_train = np.zeros((batch_repeat,Nx,Nt))
    y_train = np.zeros((batch_repeat,Nx,Nt))
    for batch in range(batch_repeat-1):
        x_train[batch,:,:] = x_snapshot
        y_train[batch,:,:] = y_snapshot
    
    return x_train, y_train, x, t