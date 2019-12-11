#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:05:41 2019

@author: morita
"""
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
import sys

from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib import rc
rc('text', usetex=True)
plt.rcParams['font.family'] = 'Times New Roman'

# %%

def load_grid(path2run,casename,Nx,Ny,Nz):
    
    ndata = Nx*Ny*Nz
    
    # cell centres
    b = 'sed \'1,22d\' %s/%s/0/C | ' %(path2run,casename)
    c = 'head -%s | sed -e \'s/(//g\' | sed -e \'s/)//g\' > ./Cdata' % (ndata)
    a = b+c
    try:
        print(a)
        subprocess.check_call(a, shell=True)
    except:
        print("Error: coundn't load 0/C")
        sys.exit(1)
    
    grid = np.loadtxt('./Cdata')

    xc = grid[:Nx,0]
    yc = grid[:,1]
    yc = yc.reshape([Ny,Nx])
    # grid data structure
    # Nx*Ny*Nz
        
    # points
    b='sed \'1,20d\' %s/%s/constant/polyMesh/points | ' %(path2run,casename)
    c='head -%s | sed -e \'s/(//g\' | sed -e \'s/)//g\' > ./pointsdata'\
        % (int((Nx+1)*(Ny+1)*(Nz+1)))
    a= b+c
    print(a)
    try:
        subprocess.check_call(a, shell=True)
    except:
        print("Error: coundn't load polyMesh/points")
        sys.exit(1)
    
    grid=np.loadtxt('./pointsdata')
    x = grid[0:Nx+1,0]
    y = np.concatenate([grid[:int((Nx+1)*(Ny/2+1)),1],
        grid[int((Nx+1)*(Ny/2+1)*(Nz+1)):int((Nx+1)*(Ny/2+1)*(Nz+1)+(Nx+1)*(Ny/2)),1]])
    y = y.reshape([Ny+1,Nx+1])
    
    print("delete temporary files")
    subprocess.check_call('rm ./Cdata', shell=True)
    subprocess.check_call('rm ./pointsdata', shell=True)
    
    return xc, yc, x, y

# %% return law data

def load_data(path2run,casename,Nx,Ny,Nz,t):
    
    ndata = Nx*Ny*Nz
    
    datalist = ('Ux','Uy','p','nut','k','omega')
    # make output file name
    num = len(datalist)
    data = np.zeros((num,ndata)) # dataname, maindata
    
    tail = np.array(["data%d" % t ]*num)
    tmp1 = np.core.defchararray.add(datalist,tail)
    for i in range(num):
        b = 'sed \'1,22d\' %s/%s/%d/%s | ' % (path2run,casename,t,datalist[i]) # delete 1-22 rows
        c = 'head -%d > ./%s' % (ndata, tmp1[i])
        a = b+c
        print(a)
        try:
            subprocess.check_call(a, shell=True)
            data[i,:] = np.loadtxt('./%s' % tmp1[i])
            subprocess.check_call('rm ./%s' % tmp1[i], shell=True)
        except:
            print("Error: Executed ReconstructPar?")
            sys.exit(1)

    Ux = data[0,:]
    Uy = data[1,:]
    p = data[2,:]
    nut = data[3,:]
    k = data[4,:]
    omega = data[5,:]
    
    Ux = Ux.reshape([Ny,Nx])
    Uy = Uy.reshape([Ny,Nx])
    p = p.reshape([Ny,Nx])
    nut = nut.reshape([Ny,Nx])
    k = k.reshape([Ny,Nx])
    omega = omega.reshape([Ny,Nx])
    # make output file name
    datalist = 'wallShearStress'
    
    tmp1 = '%sdata%d' % (datalist,t)
    b = 'sed \'1,29d\' %s/%s/%d/%s | ' % (path2run,casename,t,datalist) # delete 1-29 rows
    c = 'head -%d | sed -e \'s/(//g\' | sed -e \'s/)//g\' > ./%s' % (Nx, tmp1)
    a = b+c
    print(a)
    try:
        subprocess.check_call(a, shell=True)
    except:
        print("Error: coundn't load %d/%s" % (t,datalist))
        sys.exit(1)
    wallShearStress = np.loadtxt('./%s' % tmp1)
    tau_w = np.abs(wallShearStress[:,0]) # T_12
    subprocess.check_call('rm ./%s' % tmp1, shell=True)
    
    return Ux, Uy, p, nut, k, omega, tau_w

# %% ##################### CHECK delta99 calc. #######################
def bl_calc(Nx,Ny,Nz,xc,yc,x,y,U_infty,nu,U,V,p,nut,k,omega,tau_w):
    # nondimentionalize
    # wall unit
    u_tau = np.sqrt(tau_w)
    
    xc_p = xc*u_tau/nu
#    yc_p = np.outer(yc,u_tau)/nu
#    yp = np.outer(y,u_tau)/nu
    yc_p = np.zeros([Ny,Nx])
    yp = np.zeros([Ny+1,Nx+1])
    for i in range(Nx):
        yc_p[:,i] = yc[:,i]*u_tau[i]/nu
    
    for i in range(Nx+1):
        if i == 0:
            yp[:,0] = y[:,0]*u_tau[0]/nu
        elif i == Nx:
            yp[:,Nx] = y[:,Nx]*u_tau[Nx-1]/nu
        else:
            yp[:,i] = y[:,i]*0.5*(u_tau[i-1]+u_tau[i])/nu
#        print(yp[:,i])
    Up = np.zeros_like(U)
    kp = np.zeros_like(k)
    for i in range(Nx):
        Up[:,i] = U[:,i]/u_tau[i]
        kp[:,i] = k[:,i]/u_tau[i]**2
    
    # normal unit
    # delta99
    delta99 = np.zeros(Nx)
    U_max = np.max(U/U_infty,axis=0)
    for i in range(Nx):
        U_tmp = U[:,i]/U_max[i]
#        U_tmp = U[:,i]/U_infty # ONLY FOR DNS inflow !!!! (because of U_max fluc.)
        thre = 0.99 # initial threshhold, can be modified
        counter = 0 # avoid infinity loop
        while True:
            while True:
                for j in range(Ny):
                    if U_tmp[j] > thre:
                        index = j
                        break
                try:
#                    print('counter =',counter,',i =', i,',j =', j)
#                    print('U_tmp[index-1] =',U_tmp[index-1])
                    f1 = interpolate.interp1d(U_tmp[:index+1].reshape(-1),\
                                              yc[:index+1,i], kind="quadratic")
                    break
                except: # if U_tmp[:index] is not monotonous increase
                    thre = thre-0.001
                    counter += 1
                    if counter > 1000:
                        print('too many loop in delta99 1')
                        print('i =', i,',j =', j)
                        sys.exit(1)
            try:
                delta99[i] = f1(0.99)
#                print('counter =',counter,',i =', i,',j =', j)
#                print('U_tmp[index] =',U_tmp[index])
                break
            except: # if 0.99 is the out of the range
                thre = thre+0.0003
                counter += 1
                if counter > 1000:
                        print('too many loop in delta99 2')
                        print('i =', i,',j =', j)
                        sys.exit(1)
    
    # integral quantities
    deltaStar = np.zeros(Nx)
    theta = np.zeros(Nx)
    for i in range(Nx):
        index = np.where(yc[:,i] < delta99[i])
        U_U_infty = np.append(U[index,i]/U_max[i],0.99) # add last point
#        U_U_infty = np.append(U[index,i],0.99) # add last point
        deltaStar[i] = integrate.simps(1-U_U_infty,\
                         np.append(yc[index,i],delta99[i]))
        theta[i] = integrate.simps(U_U_infty*(1-U_U_infty),\
                     np.append(yc[index,i],delta99[i]))
    
    Re_deltaStar = U_infty*deltaStar/nu
    Re_theta = U_infty*theta/nu
    H12 = deltaStar/theta
    
    beta = np.zeros(Nx-1)
    dpdx = np.diff((p[int(Ny/2-1),:]+p[int(Ny/2),:])/2)/np.diff(xc) # middle of the channel
    for i in range(int(Nx-1)):
        beta[i] = np.mean(deltaStar[i]+deltaStar[i+1])/np.mean(tau_w[i]+tau_w[i+1])*dpdx[i]
    
    # msc.
    Cf = tau_w/(1/2*U_infty**2) # incompressible
    Re_tau = u_tau*delta99/nu
    
    return u_tau,xc_p,yc_p,yp,Up,kp,delta99,U_max,deltaStar,theta,Re_deltaStar,\
            Re_theta,H12,beta,Cf,Re_tau
    
    #  DON'T RUN TWICE!
#    U=U/U_infty
#    V=V/U_infty
#    p=p/U_infty**2
#    k=k/U_infty**2
#    nut=nut/nu
#    
#    for i in range(Nx):
#        omega[:,i]=omega[:,i]*delta99[i]/U_infty

#    return u_tau,xc_p,yc_p,yp,Up,kp,delta99,U_max,deltaStar,theta,Re_deltaStar,\
#            Re_theta,H12,beta,Cf,Re_tau,U,V,p,k,omega,nut
            
# %% misc
# ref https://qiita.com/icchi_h/items/fc0df3abb02b51f81657
def getNearestValueIndex(list, num):
    idx = np.abs(np.asarray(list) - num).argmin()
    return idx

# np vector, np vector, np vector(2D)
def distance(x,y,a): # line y=f(x) and [a0,a1]
    if np.size(x)!=np.size(y) or np.size(a)!=2:
        print('Error in distance func.: check input')
    f = interpolate.interp1d(x, y, kind="quadratic")
#    for i,val in enumerate(1,x): # i-1
#        if val>a[0]:
#            x1=x[i-1]
#            y1=y[i-1]
#            x2=val
#            y2=y[i]
#            break
    distance = np.abs(a[1]-f(a[0]))
    return distance

def getNu(path2run,casename):
    with open("%s/%s/constant/transportProperties" % (path2run,casename)) as f:
        s_line = f.readlines()
    nu = float(s_line[19][33:-2]) # NEED TO BE TUNED
    print("################### nu = %g ##################" % nu)
    return nu

    # %% some figures
    
def U_contour(casename,xc_delta,yc_delta,U,save,formatname): # give only one case
    plt.rcParams['font.size'] = 15
    figname='U_contour'
    Ny=np.shape(yc_delta)[0]
#    X, Y = np.meshgrid(xc,yc) # NY*NX
    X=np.outer(np.ones(Ny),xc_delta)
    Y=yc_delta
    
    fig = plt.figure(figsize=(12,3))
    ax = fig.add_subplot(111)
    im=ax.pcolormesh(X,Y,U, cmap='jet',vmin=0, vmax=1)
#    ax.set_aspect(3)
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
#    cax = divider.new_horizontal(size="2%", pad=0.05)
#    fig.add_axes(cax)
    clb=plt.colorbar(im,cax=cax)
    clb.set_label(r"$U/U_e^{\rm in}$")
    ax.set_xlim(0,1000)
    ax.set_ylim(0,60)
    ax.set_xlabel(r"$x/ \delta_{99}^{\rm in}$")
    ax.set_ylabel(r"$y/ \delta_{99}^{\rm in}$")
#    save_figure('U_contour','FLOW_catalog','png')
#    plt.figure()
#    plt.pcolormesh(X,Y,U, cmap='jet')
#    plt.grid()
#    pp=plt.colorbar(orientation="horizontal")
#    pp.set_label(r"$U$", fontsize=10)
#    #plt.gca().set_aspect('equal', adjustable='box')
#    plt.xlabel(r'$x$', fontsize=17)
#    plt.ylabel(r'$y$', fontsize=17)
    if save>0:
        save_figure(figname,casename,formatname)

def x_Ux(casename,Ny,x_c,Ux,save,formatname):
    figname='x_Ux'
    plt.figure()
    plt.plot(x_c,(Ux[int(Ny/2-1),:]+Ux[int(Ny/2),:])/2,color='r')
    plt.grid()
    #plt.ylim(1.1, 1.13)
    plt.xlabel(r'$x$', fontsize=17)
    plt.ylabel(r'$U_x$', fontsize=17)
    plt.title(r'$y=1$', fontsize=17)
    if save>0:
        save_figure(figname,casename,formatname)

def x_U(savename,casename,Ny,x_c,Ux,formatname,save,label):
    # up to 3 data so far
    figname='x-U'
    color=['r','b','g']
    
    n=np.size(Ny)
    plt.figure()
    for i in range(n):
        plt.plot(x_c[i,:],(Ux[i,:,:][int(Ny[i]/2-1),:]+Ux[i,:,:][int(Ny[i]/2),:])/2,\
                 color=color[i],label=label[i])
    plt.grid()
    plt.xlabel(r'$x$', fontsize=17)
    plt.ylabel(r'$U$', fontsize=17)
    plt.title(r'$y=1$', fontsize=17)
    plt.legend()
    if save>0:
        save_figure(figname,casename,formatname)

# %%
def save_figure(figname,casename,formatname):
    if 'png' in formatname:
        print('save ./figures/%s_%s.png' % (figname,casename))
        plt.savefig('./figures/%s_%s.png' % (figname,casename),bbox_inches="tight",\
                    dpi=300, pad_inches=0.05)
    if 'eps' in formatname:
        print('save ./figures/%s_%s.eps' % (figname,casename))
        plt.savefig('./figures/eps/%s_%s.eps' % (figname,casename),bbox_inches="tight",\
                     transparent=False,dpi=300, pad_inches=0.05)
    if 'pdf' in formatname:
        print('save ./figures/%s_%s.pdf' % (figname,casename))
        plt.savefig('./figures/pdf/%s_%s.pdf' % (figname,casename),bbox_inches="tight",\
                    dpi=300, pad_inches=0.05)
