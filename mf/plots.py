import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inter
import matplotlib as mpl
# from fenics import *

def connect_squares2D_gui(x,y,N,idx,connect, plotInfo, ax):
    cx, cy = ([0]*N**2 for i in range(2))
    for ii in range(0,connect.shape[0]):
        cx[ii] = np.mean(x[0,idx[connect[ii,:]]])
        cy[ii] = np.mean(y[0,idx[connect[ii,:]]])
        
        if plotInfo in 'yes':
            ax.plot( [x[0, idx[connect[ii,0]] ], x[0, idx[connect[ii,1]] ] ], [y[0, idx[connect[ii,0]] ], y[0, idx[connect[ii,1]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            ax.plot( [x[0, idx[connect[ii,1]] ], x[0, idx[connect[ii,2]] ] ], [y[0, idx[connect[ii,1]] ], y[0, idx[connect[ii,2]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            ax.plot( [x[0, idx[connect[ii,2]] ], x[0, idx[connect[ii,3]] ] ], [y[0, idx[connect[ii,2]] ], y[0, idx[connect[ii,3]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            ax.plot( [x[0, idx[connect[ii,3]] ], x[0, idx[connect[ii,0]] ] ], [y[0, idx[connect[ii,3]] ], y[0, idx[connect[ii,0]] ] ], ':', color=[.6,.6,.6], lw=0.5)
        
        cs = np.c_[np.array(cx),np.array(cy)]
    return cs

def create_plate2D_gui(Lx,Ly,N,plot_info, ax):
    x = np.linspace(0,Lx,N+1)    # 0 : dx : Lx ;              
    y = np.linspace(0,Ly,N+1)    # 0 : dy : Ly ;
    
    X,Y = np.meshgrid(x,y)
    xx = X.reshape(1, (N+1)**2)
    yy = Y.reshape(1, (N+1)**2)
    
    idx = np.array(np.r_[:xx.shape[1]],dtype='int')
    connect =np.zeros([4,N**2], dtype='int')
    
    for ii in range(N):
        for jj in range(4):
            connect[0,N*ii:N*(ii+1)] =  np.r_[(ii+0) * N+1 + ii: (ii+1) * N+1 + ii]
            connect[1,N*ii:N*(ii+1)] =  np.r_[(ii+0) * N+2 + ii: (ii+1) * N+2 + ii]
            connect[2,N*ii:N*(ii+1)] =  np.r_[(ii+1) * N+3 + ii: (ii+2) * N+3 + ii]
            connect[3,N*ii:N*(ii+1)] =  np.r_[(ii+1) * N+2 + ii: (ii+2) * N+2 + ii]
    
    connect = connect.T - 1
    cs = connect_squares2D_gui(xx, yy, N,idx, connect, plot_info, ax)
    
    return cs, x, y

def connect_squares2D(x,y,N,idx,connect, plotInfo):
    cx, cy = ([0]*N**2 for i in range(2))
    for ii in range(0,connect.shape[0]):
        cx[ii] = np.mean(x[0,idx[connect[ii,:]]])
        cy[ii] = np.mean(y[0,idx[connect[ii,:]]])
        
        if plotInfo in 'yes':
            plt.plot( [x[0, idx[connect[ii,0]] ], x[0, idx[connect[ii,1]] ] ], [y[0, idx[connect[ii,0]] ], y[0, idx[connect[ii,1]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            plt.plot( [x[0, idx[connect[ii,1]] ], x[0, idx[connect[ii,2]] ] ], [y[0, idx[connect[ii,1]] ], y[0, idx[connect[ii,2]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            plt.plot( [x[0, idx[connect[ii,2]] ], x[0, idx[connect[ii,3]] ] ], [y[0, idx[connect[ii,2]] ], y[0, idx[connect[ii,3]] ] ], ':', color=[.6,.6,.6], lw=0.5)
            plt.plot( [x[0, idx[connect[ii,3]] ], x[0, idx[connect[ii,0]] ] ], [y[0, idx[connect[ii,3]] ], y[0, idx[connect[ii,0]] ] ], ':', color=[.6,.6,.6], lw=0.5)
        
    cs = np.c_[np.array(cx),np.array(cy)]
    return cs

def create_plate2D(Lx,Ly,N,plot_info):
    x = np.linspace(0,Lx,N+1)    # 0 : dx : Lx ;              
    y = np.linspace(0,Ly,N+1)    # 0 : dy : Ly ;
    
    X,Y = np.meshgrid(x,y)
    xx = X.reshape(1, (N+1)**2)
    yy = Y.reshape(1, (N+1)**2)
    
    idx = np.array(np.r_[:xx.shape[1]],dtype='int')
    connect =np.zeros([4,N**2], dtype='int')
    
    for ii in range(N):
        for jj in range(4):
            connect[0,N*ii:N*(ii+1)] =  np.r_[(ii+0) * N+1 + ii: (ii+1) * N+1 + ii]
            connect[1,N*ii:N*(ii+1)] =  np.r_[(ii+0) * N+2 + ii: (ii+1) * N+2 + ii]
            connect[2,N*ii:N*(ii+1)] =  np.r_[(ii+1) * N+3 + ii: (ii+2) * N+3 + ii]
            connect[3,N*ii:N*(ii+1)] =  np.r_[(ii+1) * N+2 + ii: (ii+2) * N+2 + ii]
    
    connect = connect.T - 1
    cs = connect_squares2D(xx, yy, N,idx, connect, plot_info)
    
    return cs, x, y