#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 15:11:58 2019

@author: morita
"""
import numpy as np
import sys

import postProcess_func

# %% logging
import logging
# # create logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
# logger.setLevel(logging.INFO)

# %% global variables
PATH2RUN ='../..'
caseDict = \
    {"prec":{"Uinf":1,"delta99_in":0.05,"Nx":5000,"Ny":218,"Nz":1,"t":198},\
     "prec_xssscoarse":{"Uinf":1,"delta99_in":0.05,"Nx":500,"Ny":218,"Nz":1,"t":198},\
     "prec_xsscoarse":{"Uinf":1,"delta99_in":0.05,"Nx":625,"Ny":218,"Nz":1,"t":198},\
     "prec_xscoarse":{"Uinf":1,"delta99_in":0.05,"Nx":1250,"Ny":218,"Nz":1,"t":198},\
     "prec_xcoarse":{"Uinf":1,"delta99_in":0.05,"Nx":2500,"Ny":218,"Nz":1,"t":198},\
     "prec_xfine":{"Uinf":1,"delta99_in":0.05,"Nx":10000,"Ny":218,"Nz":1,"t":198},\
     "prec_ycoarse":{"Uinf":1,"delta99_in":0.05,"Nx":5000,"Ny":108,"Nz":1,"t":198},\
     "prec_yfine":{"Uinf":1,"delta99_in":0.05,"Nx":5000,"Ny":440,"Nz":1,"t":198},\
     "main_1wall":{"Uinf":1,"delta99_in":0.05,"Nx":5000,"Ny":238,"Nz":1,"t":198},\
     "GPtest":{"Uinf":1,"delta99_in":0.05,"Nx":2500,"Ny":238,"Nz":1,"t":120},\
     "GPtest_xcoarse":{"Uinf":1,"delta99_in":0.05,"Nx":1250,"Ny":238,"Nz":1,"t":120},\
     }

# %% functions
def get_bl_ref(year):
    if year == 2009:
        # from DNS <https://www.mech.kth.se/~pschlatt/DATA/vel_0670.prof>
        # Turbulent boundary layers up to Re_theta=2500 studied through simulation and experiment
        Re_theta = np.array([670.846, 1000.490,1412.843,2000.404,2400.206])
        Re_deltaStar = np.array([1010.400,1467.762,2032.913, 2837.624,3384.808])
        Re_tau = np.array([254.7530,365.4094,500.7734,683.5205, 802.2587])
        H12 = np.array([1.506158, 1.467043,1.438881,1.418526,1.410216])
        cf = np.array([0.004740071,0.004258940,0.003893054  ,0.003556352  , 0.003388410])
        
    elif year == 2010:
        #  from DNS <https://www.mech.kth.se/~pschlatt/DATA/vel_0670_dns.prof>
        Re_theta = \
        np.array([677.452,1006.534,1420.960,2000.703,2536.827,3031.843,3273.634,3626.304,3969.154,4061.378])
        Re_deltaStar = \
        np.array([998.106,1459.397,2030.877,2827.938,3563.325,4237.594,4567.562,5044.407,5508.503,5633.318])
        Re_tau = \
        np.array([252.2550,359.3794,492.2115,671.1240,830.0115,974.1849,1043.4272,1145.1699,1244.7742,1271.5350])
        H12 = \
        np.array([1.473324,1.449924,1.429229,1.413472,1.404639,1.397696,1.395257,1.391060,1.387828,1.387046])
        cf = \
        np.array([0.004777047,0.004264437,0.003884512,0.003539148,0.003327997,0.003189402,\
               0.003123259,0.003055599,0.002987022,0.002970989])
        
    elif year == 2014:
        # Eitel-Amor, Ramis and Schlatter, Int. J. Heat Fluid Flow 47, 57-69, 2014
        Re_theta = \
        np.array([1375.753,2240.616,3014.185,3739.158,4427.882,5095.319,5746.112,\
                  6381.012,7000.797,7603.115,8183.195])
        Re_deltaStar = \
        np.array([1956.086,3138.983,4194.470,5177.661,6109.742,7002.279,7865.991,\
                  8705.101,9519.846,10307.759,11065.409])
        Re_tau = \
        np.array([457.8235,725.4953,956.9811,1169.2280,1367.3586,1561.0620,1750.5198,\
                  1937.3113,2118.0861,2299.2119,2478.9901])
        H12 = \
        np.array([1.421829,1.400946,1.391577,1.384713,1.379834,1.374257,1.368924,\
                  1.364219,1.359823,1.355728,1.352211])
        cf = \
        np.array([0.004002402,0.003513976,0.003252717,0.003084523,0.002957985,\
                  0.002871772,0.002803321,0.002746413,0.002695571,0.002653308,0.002623404])
        
    else:
        logger.error("Given year should be either 2009, 2010 or 2014")
        sys.exit(1)
    
    return Re_theta,Re_deltaStar,Re_tau,H12,cf
#    u_tau = np.sqrt(1/2*U_infty**2*cf)
#    delta99 = Re_tau*nu/u_tau
#    deltaStar = Re_deltaStar*nu/U_infty
#    theta = Re_theta*nu/U_infty
#    Re_delta=U_infty*delta99/nu
#    
#    return Re_theta,Re_deltaStar,Re_tau,H12,cf,u_tau,\
#        delta99,deltaStar,theta,Re_delta

# from 2010(DNS) or 2014(LES) data
def get_refProfile(path2file,LESflag=False): # default DNS
    if LESflag:
        skiprows=12
    else:
        skiprows=14
        
    try:
        tmp = np.loadtxt(path2file, skiprows=skiprows, unpack=True)
        # y_delta99, yp, Up, urms_p, vrms_p, wrms_p\
        #     = np.loadtxt(path2file, skiprows=skiprows, unpack=True)[:6]
        # Vp = np.loadtxt(path2file, skiprows=skiprows)
    except:
        logger.error("Cannot read %s" % path2file)
        sys.exit(1)
    
    y_delta99, yp, Up, urms_p, vrms_p, wrms_p = tmp[:6]
    if LESflag:
        Vp = tmp[:,-1]
    else:
        Vp = tmp[:,13]
    
    return y_delta99, yp, Up, Vp, urms_p, vrms_p, wrms_p

# def get_refProfileAll(LESflag=0):
#     if LESflag:
#         Re_thetaList = np.array(['01000','02000','03000','04000','05000','06000',\
#                          '07000','08000','09000','10000','11000'])
#     else:
#         Re_thetaList = np.array(['0670','1000','1410','2000','2540','3030','3270',\
#                              '3630','3970','4060'])
    
#     y_delta = np.empty(0)
#     yp = np.empty(0)
#     Up = np.empty(0)
    
#     for i in Re_thetaList:
#         if LESflag:
#             path2file="./bl_data/LES/RE8000/vel_%s_.prof" % i
#         else:
#             path2file="./bl_data/vel_%s_dns.prof" % i
#         tmp,tmp1,tmp2 = get_refProfile(path2file,LESflag)[:3]
#         # 1D array
#         y_delta = np.append(y_delta,tmp)
#         yp = np.append(yp,tmp1)
#         Up = np.append(Up,tmp2)
    
#     y_delta = y_delta.reshape([len(Re_thetaList),-1])
#     yp = yp.reshape([len(Re_thetaList),-1])
#     Up = Up.reshape([len(Re_thetaList),-1])
#     return y_delta, yp, Up

# %% classes
class RefData(): 
    def __init__(self,year): # get integral quantities
        if year not in [2009, 2010, 2014]:
            logger.critical("Input error")
            sys.exit(1)
        self.year = year
        self.Re_theta, self.Re_deltaStar, self.Re_tau, self.H12, self.cf\
            = get_bl_ref(self.year)
    
    def get_profile(self, Re_theta): # only for 2010 & 2014, Re_theta should be string
        if self.year != 2010 and self.year != 2014:
            logger.critical("profile available only from 2010, 2014")
            sys.exit(1)
            
        if self.year == 2014:
            LESflag = True
        else:
            LESflag = False
            
        if LESflag:
            path2file="./bl_data/LES/RE8000/vel_%s_.prof" % Re_theta
        else:
            path2file="./bl_data/vel_%s_dns.prof" % Re_theta
            
        self.y_delta, self.yp, self.Up, self.Vp, self.urms_p, self.vrms_p, self.wrms_p\
            = get_refProfile(path2file, LESflag)
    
    # def get_profileAll(self): # only for 2010 & 2014
    #     if self.year != 2010 and self.year != 2014:
    #         logger.critical("profile available only from 2010, 2014")
    #         sys.exit(1)
    #     if self.year == 2014:
    #         LESflag = 1
    #     else:
    #         LESflag = 0
    #     self.y_delta, self.yp, self.Up = get_refProfileAll(LESflag)
    
    def set_label(self, labelName):
        self.label = labelName
    
    def set_marker(self, marker): # ex. "C1+"
        self.marker = marker
    
    def set_line(self, line): # ex. "C1:"
        self.line = line


class OFresult():
    def __init__(self, caseName, caseConfig, path2run=PATH2RUN, autoLoad=True):
        self.path2run = path2run
        self.caseName = caseName
        self.dict = caseConfig
        self.dict["nu"] = postProcess_func.getNu(self.path2run, self.caseName)
        if autoLoad==True:
            logger.info("start auto load")
            self.get_grid()
            self.get_result()
            self.calc()
        logger.debug("class: OFresult is constructed")
    
    # for plot
    def set_label(self, labelName):
        self.label = labelName
    
    def set_marker(self, marker): # ex. "C1+"
        self.marker = marker
    
    def set_line(self, line): # ex. "b:"
        self.line = line
    
    # method can be used to overwrite
    def get_grid(self):
        self.xc, self.yc, self.x, self.y\
            = postProcess_func.load_grid(self.path2run, self.caseName, self.dict["Nx"],\
                                         self.dict["Ny"], self.dict["Nz"])
                
    def get_result(self):
        self.U, self.V, self.p, self.nut, self.k, self.omega, self.tau_w\
            = postProcess_func.load_data(self.path2run, self.caseName, self.dict["Nx"],\
                                         self.dict["Ny"], self.dict["Nz"], self.dict["t"])
                
    def calc(self):
        self.u_tau, self.xc_p, self.xp, self.yc_p, self.yp, self.Up, self.kp, self.delta99,\
            self.U_max, self.deltaStar, self.theta, self.Re_deltaStar, self.Re_theta,\
            self.H12, self.beta, self.cf, self.Re_tau\
                = postProcess_func.bl_calc(self.dict["Nx"], self.dict["Ny"], self.dict["Nz"],\
                                           self.xc, self.yc, self.x, self.y, self.dict["Uinf"],\
                                           self.dict["nu"], self.U, self.V, self.p, self.nut,\
                                           self.k, self.omega, self.tau_w)
        self.xc_delta = self.xc/self.dict["delta99_in"]
        self.x_delta = self.x/self.dict["delta99_in"]
        self.yc_delta = self.yc/self.dict["delta99_in"]
        self.y_delta = self.y/self.dict["delta99_in"]
        self.Re_delta = self.dict["Uinf"]*self.delta99/self.dict["nu"]
        self.Re_x = self.dict["Uinf"]*self.xc/self.dict["nu"]
        # Re_theta corresponds to xc
        self.Re_theta_c = \
            np.array([(self.Re_theta[i]+self.Re_theta[i+1])/2 for i in range(self.dict["Nx"]-1)])
            
    def assess_grid(self):
        print(self.caseName)
        max_xp = max(np.diff(self.xp))
        mean_xp = np.mean(np.diff(self.xp))
        max_yp = max(self.yp[1,:])
        mean_yp = np.mean(self.yp[1,:])
        print("max(x+) =", max_xp)
        print("mean(x+) =", mean_xp)
        print("max(y+) =", max_yp)
        print("mean(y+) =", mean_yp)
