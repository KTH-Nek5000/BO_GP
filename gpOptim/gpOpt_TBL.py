##################################################################################
#*********** Bayesian Optimization using Gaussian Processes ****************
# Find optimal shape of the upper boundary s.t. th TBL at the lower wall has a 
#      specific value.
##################################################################################
# Saleh Rezaeiravesh, salehr@kth.se
#---------------------------------------------------------------------------------

import sys
import os
import math as mt
import matplotlib
matplotlib.use('PDF')
# import matplotlib.cm as cm
# import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
plt.rcParams["font.size"] = 20
rc('text', usetex=True)
import numpy as np
import GPy
import GPyOpt
from GPyOpt.methods import BayesianOptimization
# from GPyOpt import Design_space
# from GPyOpt.experiment_design import initial_design
from numpy.linalg import norm
from matplotlib.colors import Normalize

# %% logging
import logging
# # # create logger
# logger = logging.getLogger("gpOptim/gpOpt_TBL.py") # root logger
# if (logger.hasHandlers()):
#     logger.handlers.clear()
# logger.setLevel(logging.INFO)

# def add_handler():
#     # create console handler and set level to debug
#     ch = logging.StreamHandler()
#     ch.setLevel(logging.INFO)
#     formatter = logging.Formatter('%(name)s - %(funcName)s - %(levelname)s - %(message)s')
#     ch.setFormatter(formatter)
#     # if not logger.handlers:
#     #     logger.addHandler(ch)
#     logger.addHandler(ch)

# add_handler()

logger = logging.getLogger("Driver").getChild("gpOptim/gpOpt_TBL.py")

# %% global variables
#----------------------------------------------------------------------------
#>>>> SETTINGS & PROBLEM DEFINITION -----------------------------------------
sigma_d=0.0       #sdev of the white noise in the measured data   
whichOptim='min'  #find 'max' or 'min' of f(x)?
# tol_d=0.02        #minimum distace between two consequtive samples x to keep the code running
# tol_b=0.1        #deviation between best f(x+) in two consequtive iterations (relative error)
                  #note if err_d<tol_d and err_b<tol_b => convergence in (x_opt , f(x_opt))
tol_abs=0.1
kernelType='Matern52'  #'RBF', 'Matern52'
#admissible range of parameters
qBound=[[40,46], [40,44]] # /delta99^in
nPar = np.shape(qBound)[0] #number of parameters, p  dimension of x={x1,x2,...,xp} where y=f(x)
nGPinit=1   #minimum number of GP samples in the list to start BO-GP algorithm
            #to avoid random sampling from the parameter space: see nextGPsample()
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

#>>>>Assignments (don't touch these!)
#admissible domain in GPy format
domain=[];
for i in range(nPar):
    domain_={'name':'q'+str(i+1),'type':'continuous','domain':(qBound[i][0],qBound[i][1])}
    domain.append(domain_)

if whichOptim=='max':
    maxFlag=True
    ifac=-1.0
else:
    maxFlag=False
    ifac=1.0
    
#evaluating the true function? (no noise)
if sigma_d==0.:
    fevalFlag=True
else:
    fevalFlag=False
    
#domain bounds
BOUNDS=[]
for i in range(nPar):
    BOUNDS.append(domain[i]['domain'])

# %% INT FUNCTIONS
def read_available_GPsamples(gpInputFile,nPar):
    """
        Read the most updated list of (x,y) GP samples from gpInputFile
    """
    F1=open(gpInputFile,'r')
    ain=F1.readlines()
    ain_sep=[];
    for i in range(len(ain)):
        ain_sep.append(ain[i].split())
    iskip=2;  # no. of lines to skip from the beginning of the input file
    n=len(ain_sep)-iskip;
    xList=np.zeros((n,nPar))
    yList=np.zeros(n)
    for i in range(n):
        for j in range(nPar):
            xList[i,j]=float(ain_sep[i+iskip][j+1])
        yList[i]=float(ain_sep[i+iskip][nPar+1])
    F1.close()
    if n < 1:
        logger.warning("No available sample in %s" % gpInputFile)
    else:
        logger.info("read available samples from %s" % gpInputFile)
    return xList,yList

#//////////////////////////////////////////////////////////
def update_GPsamples(gpOutputFile,xList,yList,xNext,yNext):
    """
        Update the existing list of GP samples with the recent sample & response
    """
    F2=open(gpOutputFile,'w')
    F2.write('#List of GP samples." \n')
    tmp='#iter'+'\t'
    for i in range(len(xNext)):
        tmp=tmp+'p'+str(i+1)+'\t'
    tmp=tmp+'response\n'
    F2.write(tmp)
    nData=xList.shape[0]
    #write the datalist before observing the new sample
    for i in range(nData):
        tmpList=str(i+1)+'\t'
        for j in range(nPar):
            tmpList=tmpList+str(xList[i][j])+'\t'
        tmpList=tmpList+str(yList[i][0])+'\n'
        F2.write(tmpList)
    #write the most recent parameters and associated response
    tmpList=str((nData+1))+'\t'
    for j in range(nPar):
        tmpList=tmpList+str(xNext[j])+'\t'
    tmpList=tmpList+str(yNext)+'\n'
    F2.write(tmpList)
    F2.close()
    logger.info('**** %s is updated!' % gpOutputFile)

#//////////////////////////////////////////////////////////////
def my_convergence_plot(xList,yList,whichOptim,figDir,figName):
    """
       Plot convergence of parameter samples and associated response
    """
    xDistList=[];  #Distance between two consecutive parameter samples  
    yBestList=[];  #Best response value up to any iteration
    nData=xList.shape[0];
    for i in range(1,nData):
        xDistList.append(norm(xList[i][:]-xList[i-1][:])) # 2-norm
    for i in range(nData):
        if whichOptim=='min':
           yBestList.append(min(yList[:i+1]))
        elif whichOptim=='max':
           yBestList.append(max(yList[:i+1]))
    
    fig=plt.figure()
    plt.subplot(2,1,1)
    plt.semilogy(range(1,nData),xDistList,'-ob',lw=2)
    plt.title("Distance between 2 consecutive parameter samples",fontsize=20)
    plt.xlabel(r'$N$',fontsize=20)
    plt.ylabel(r'$\| \mbox{\boldmath{$q$}}^{(n+1)} - \mbox{\boldmath{$q$}}^{(n)} \|_2$',fontsize=22)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True)
    plt.gca().get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
    plt.subplot(2,1,2)
    plt.semilogy(range(1,nData+1),yBestList,'-or',lw=2)
    plt.title('Best Value So Far',fontsize=20)
    plt.xlabel(r'$N$',fontsize=20)
    plt.ylabel(r'${\rm min}~(\mathcal{R})$',fontsize=22)
    plt.gca().get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
    plt.tick_params(labelsize=20)
    plt.grid(True)
    fig = plt.gcf()
    DPI = fig.get_dpi()
    fig.set_size_inches(600/float(DPI),1200/float(DPI))
    plt.savefig(figDir+figName+'.pdf',bbox_inches='tight')
    logger.info('save: %s%s.pdf' % (figDir,figName))
    #plt.show()
    return xDistList,yBestList

#////////////////////////////////////////////
def test_grid(bounds1,bounds2,nTest1,nTest2):
    """
       Construct a 2D mesh (test data) with uniformly spaced points to 
       illustrate contours of the response
    """
    nTest=nTest1*nTest2
    x1Test=np.linspace(bounds1[0],bounds1[1],nTest1)
    x2Test=np.linspace(bounds2[0],bounds2[1],nTest2)
    x1TestGrid=np.zeros((nTest1,nTest2))  
    x2TestGrid=np.zeros((nTest1,nTest2))
    # yTestGrid=np.zeros((nTest1,nTest2))
    xTestArr=np.zeros((nTest,2));
    for i in range(nTest1):
        for j in range(nTest2):
            k=i+j*nTest1
            xTestArr[k,0]=x1Test[i]
            xTestArr[k,1]=x2Test[j]
            x1TestGrid[i,j]=x1Test[i]
            x2TestGrid[i,j]=x2Test[j]
    xTestArr=np.asarray(xTestArr)   #n* x p=2
    return x1TestGrid,x2TestGrid,xTestArr

#>>>>> plots
#////////////////////////////////////////////////////////////////////
def gpOpt1d_postProc(xGP,yGP,sigma_d,domain,nTest,plotOpts):
    """ 
       Postprocess the Bayesian optimization on a 1D parameter space.
       This method constructs a GPR based on the final set of GP samples taken 
           during the optimization. The hyper-parameters of the GPR are optimized. 
       NOTE: only for nPar==1 case
    """
    #GP kernel
    if plotOpts['kernelType']=='RBF':
       K=GPy.kern.RBF(input_dim=1,lengthscale=1.0,variance=1.0)
    elif plotOpts['kernelType']=='Matern52':
       K=GPy.kern.Matern52(input_dim=1,lengthscale=1.0,variance=1.0)
    #define the GPR
    gprFinal=GPy.models.GPRegression(xGP,yGP,kernel=K,noise_var=sigma_d**2.)
    gprFinal.constrain_positive()  #make all parameters positive

    #if you want to get exactly the same plot as "GPyOpt.plot_acquisition()", 
    # you need to let gaussicna_noise.variance to be optimized (unfixed) 
    gprFinal.Gaussian_noise.variance.fix()   #sigma_d = fixed
    gprFinal.optimize('bfgs', max_iters=200)   #optimization of hyperparameters

    #make predictions at test points
    domain_=domain[0]['domain']
    xTest_=np.linspace(domain_[0],domain_[1],nTest)
    xTest=xTest_.reshape(nTest,1)
    [meanPred,covarPred]=gprFinal.predict(xTest,full_cov=True)
    gpyPlotter_1D(meanPred[:,0],covarPred,xGP,yGP,xTest,plotOpts)

#//////////////////////////////////////////////////////////
def gpOpt2d_postProc(nPar,xGP,yGP,sigma_d,bounds,plotOpts):
    """ 
       Postprocess the Bayesian optimization on a 2D parameter space.
       This method constructs a GPR based on the final set of GP samples taken 
           during the optimization. The hyper-parameters of the GPR are optimized. 
       NOTE: for Now nPar should be either 2 or 4. 
    """
    #>>> 0. Assign the ID od mutual parameters
    parID=[]
    if nPar==2:
        parID.append([0,1])
    elif nPar==3:
        parID.append([0,1])
        parID.append([0,2])
        parID.append([1,2])
        loc=[0,2,3]
    elif nPar==4:
       parID.append([0,1])
       parID.append([0,2])
       parID.append([1,2])
       parID.append([0,3])
       parID.append([1,3])
       parID.append([2,3])
       loc=[0,3,4,6,7,8]  #location in subplot
       #
       # parID.append([1,3])
       # parID.append([0,2])
       # parID.append([1,2])
       # parID.append([3,2])
       # parID.append([1,0])
       # parID.append([3,0])
       # loc=[7,4,1,5,2,3]  #location in subplot
    else:
        logger.error("nPar should be 2, 3 or 4: given %d" % nPar)
    
    # plt.figure()
    for i in range(len(parID)):  #param-pair loop
        I=parID[i][0]   #ID of param 1 in the pair
        J=parID[i][1]   #ID of param 2

        #>>> 1. Construct the GPR for each 2 parameters
        #assign GP kernel
        if plotOpts['kernelType']=='RBF':
          K=GPy.kern.RBF(input_dim=2,lengthscale=1.0,variance=1.0)
        elif plotOpts['kernelType']=='Matern52':
          K=GPy.kern.Matern52(input_dim=2,lengthscale=1.0,variance=1.0)
        #define the GPR
        xGP_=[]
        for j in range(len(yGP)):
            xGP_.append([xGP[j][I] , xGP[j][J]])
        xGP_=np.asarray(xGP_)
        gprFinal=GPy.models.GPRegression(xGP_,yGP,kernel=K,noise_var=sigma_d**2.)
        gprFinal.constrain_positive()  #make all parameters positive
        # if you want to get exactly the same plot as "GPyOpt.plot_acquisition()",
        # you need to let gaussicna_noise.variance to be optimized (unfixed) 
        gprFinal.Gaussian_noise.variance.fix()     #sigma_d = fixed
        gprFinal.optimize('bfgs', max_iters=200)   #optimization of hyperparameters
        logger.info('------------------------------------------')
        logger.info('Final GPR with optimal hyper-parameters:')
        logger.info('--------Model paramaters (%d, %d) -----------' % (I+1, J+1))
        logger.info(str(gprFinal))

        #>>> 2. Generate a mesh of test parameters in a 2d-plane
        bound1=[bounds[I][0],bounds[I][1]]
        bound2=[bounds[J][0],bounds[J][1]]
        [x1_grid,x2_grid,x_grid]=test_grid(bound1,bound2,40,40)

        #>>> 3. Make predictions at the test samples
        [meanPred,covarPred]=gprFinal.predict(x_grid)

        #>>> 4. Plot response surface predicted by GPR at test mesh
        #plot the GPR      
        if nPar==2:
            plt.subplot(1,1,1)
            figSize=500
        elif nPar==3:
            plt.subplot(2,2,loc[i])
            figSize=1500
        elif nPar==4:
            plt.subplot(3,3,loc[i])
            figSize=1500
        gpyPlotter_2Dc(meanPred,covarPred,xGP_,yGP,x1_grid,x2_grid,I,J,plotOpts)

    #save fig
    if 'figDir' in plotOpts.keys():
       figDir=plotOpts['figDir']
       if not os.path.exists(figDir):
          os.makedirs(figDir)
    if 'figName' in plotOpts.keys():
       figName=plotOpts['figName']
       figSave=figDir+figName
    fig = plt.gcf()
    DPI = fig.get_dpi()
    fig.set_size_inches(figSize/float(DPI),figSize/float(DPI))
    plt.savefig(figSave+'.pdf',bbox_inches='tight')
    logger.info("save: %s.pdf" % figSave)

#//////////////////////////////////////////////////////////////////
def gpyPlotter_1D(meanPred,covarPred,xGP,yGP,xTest_,plotOpts):
    """ 
       plot training data (GP samples), mean prediction, 95% confidence, an abitrary sample
    """
    xTest=xTest_[:,0]
    #Assigning
    # ifac=1.0
    if 'whichOptim' in plotOpts.keys():
       if plotOpts['whichOptim']=='max':   #for plotting the GPR after computing max
          meanPred=-meanPred;
          yTrain=-yGP
          # ifac=-1.0
    n=len(meanPred);
    #confidence interval
    confidPred=[];
    for i in range(n):
        confidPred.append(1.96*mt.sqrt(covarPred[i,i]))  #95% confidence interval
    #arbitrary sample from the prediction distribution
    if 'arbitSample' in plotOpts.keys(): arbitSamp=plotOpts['arbitSample']
    if (arbitSamp=='yes' and covarPred.ndim == 2 and covarPred.shape[0]==covarPred.shape[1]):
       ySample1=np.random.multivariate_normal(meanPred,covarPred)
    else:  #if single test data is used
       ySample1=[]
    #plot
    plt.figure(figsize=(20,8));
    ax=plt.gca();
    ax.fill_between(xTest,meanPred+confidPred,meanPred-confidPred,color='powderblue',alpha=0.5)
    plt.semilogy(xTest,meanPred+confidPred,'-',color='b',linewidth=1,alpha=0.2)#,label=r'$95\%$ confidence')
    plt.semilogy(xTest,meanPred-confidPred,'-',color='b',linewidth=1,alpha=0.2)
    plt.semilogy(xGP[:-1],yGP[:-1],'o r',markersize=7,label='Training Data')
    plt.semilogy(xGP[-1],yGP[-1],'o g',markersize=10)#,label='Training Data')
    plt.plot(xTest,meanPred,'-b',linewidth=2,label='Predicted Mean')
    if (len(ySample1)>0):
       plt.plot(xTest,ySample1,'--r',linewidth=2,label='Arbitrary Sample',dashes=[3,2])

    xmin=qBound[0][0]
    xmax=qBound[0][1]
    plt.xlabel(r'$q$',fontsize=22)
    plt.ylabel(r'$\mathcal{R}$',fontsize=22)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(loc='upper left',fontsize=20);
    plt.xlim(xmin,xmax)
    if 'ylim' in plotOpts.keys():
       ylim=plotOpts['ylim']
       plt.ylim((ylim[0],ylim[1]))
    plt.title(r"$N_i = %d,~ \min(\mathcal{R}) = %f$" % (np.size(yGP),xGP[-1],np.min(yGP)),fontsize=20)
    plt.grid(True);
    #save fig
    if 'figDir' in plotOpts.keys():
       figDir=plotOpts['figDir']
       if not os.path.exists(figDir):
          os.makedirs(figDir)
    if 'figName' in plotOpts.keys():
       figName=plotOpts['figName']
       figSave=figDir+figName
    fig = plt.gcf()
    DPI = fig.get_dpi()
    fig.set_size_inches(1000/float(DPI),500/float(DPI))
    if 'figSave' in locals():
       plt.savefig(figDir+figName+'.pdf',bbox_inches='tight')
       logger.info("save figure as %s%s.pdf" % (figDir,figName))
    #plt.show()

#///////////////////////////
def gpyPlotter_2Dc(meanPred,covarPred,x,y,x1TestGrid,x2TestGrid,I,J,plotOpts):
    """ 
    Plot 2D contour lines generated by GPR + the surface of the error between
            true function and the perdictions by GPR
       (meanPred,covarPred): mean and covariance predicted by GPR at grid test points
       (x,y): GP samples and associated model function values
       (x1TestGrid,x2TestGrid): coordinates of the test grid points
       yTestGrid: values of the true function at xTestGrid
       yTestGrid: GPR prediction for the function value at (x1TestGrid,x2TestGrid)
    """
    #Assigning
    # ifac=1.0
    if 'whichOptim' in plotOpts.keys():
       if plotOpts['whichOptim']=='max':   #for plotting the GPR after computing max
          #meanPred=-meanPred;
          # ifac=-1.0
           pass
    nTest1=x1TestGrid.shape[0]
    nTest2=x1TestGrid.shape[1]
    nTest=nTest1*nTest2
    #Confidence interval
    confidPred=[];
    for i in range(nTest):
        if covarPred[i]<0:
           logger.warning('negative predicted variance = %g at i = %d is replaced by 0.0'\
                          %(covarPred[i,i], i))
           covarPred[i]=0.0
        confidPred.append(1.96*mt.sqrt(covarPred[i]))  #95% confidence interval

    #reshape arrays for plotting    
    meanPredGrid=np.zeros((nTest1,nTest2))  #only for plotting
    covarPredGrid=np.zeros((nTest1,nTest2))  #only for plotting
    for i in range(nTest1):
        for j in range(nTest2):
            k=i+j*nTest1
            meanPredGrid[i,j]=meanPred[k]
            covarPredGrid[i,j]=covarPred[k]

    #2D contourplot
    # ax=plt.gca()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if "Rmin" in plotOpts.keys():
        # CS=plt.contourf(x1TestGrid,x2TestGrid,meanPredGrid,40,cmap="jet",vmin=plotOpts["Rmin"],vmax=plotOpts["Rmax"])#,label=r'$95\%$ confidence')
        CS=ax.contourf(x1TestGrid,x2TestGrid,meanPredGrid, \
                       levels=np.linspace(plotOpts["Rmin"],plotOpts["Rmax"],41), \
                           cmap="jet",norm=Normalize(vmin=plotOpts["Rmin"], vmax=plotOpts["Rmax"]),extend='both')
        # CS.set_clim(vmin=plotOpts["Rmin"], vmax=plotOpts["Rmax"])
        clb = fig.colorbar(CS)
        # print(plotOpts["Rmax"])
        
    else:
        CS=plt.contourf(x1TestGrid,x2TestGrid,meanPredGrid,40,cmap="jet")#,label=r'$95\%$ confidence')
        clb = plt.colorbar(CS)
    # plt.clabel(CS, inline=True, fontsize=13,colors='k',fmt='%0.2f',rightside_up=True,manual=False)
    clb.set_label(r"$\mathcal{R}$")
    plt.plot(x[:,0],x[:,1],'--ok',markersize=7,label='Training Data')
    plt.plot(x[len(y)-1,0],x[len(y)-1,1],'--sr',markersize=7,label='Training Data')
    # plt.title(r'Mean GPR Prediction',fontsize=18)
    plt.xlabel(r'$%s%s$' % ("q_", str(I+1)),fontsize=20)
    plt.ylabel(r'$%s%s$' % ("q_", str(J+1)),fontsize=20)
    ##contours of uncertainty 
    ##plot only if nPar==2
    #plt.subplot(1,2,2)
    #ax=plt.gca()
    #CS=plt.contour(x1TestGrid,x2TestGrid,1.96*np.sqrt(covarPredGrid),40)#,label=r'$95\%$ confidence')
    #plt.clabel(CS, inline=True, fontsize=13,colors='k',fmt='%0.2f',rightside_up=True,manual=False)
    #plt.xlabel('q'+str(I+1),fontsize=20)
    #plt.ylabel('q'+str(J+1),fontsize=20)
    #plt.title(r'$95\%$ Uncertainty',fontsize=18)
    # fig.tight_layout()
    return fig

# %% EXECUTABLE FUNCTIONS
def printSetting():
    """
    print global parameters
    """
    logger.info("\nnPar = %d\nsigma_d = %f\nwhichOptim = %s\ntol_abs = %f"\
                "\nkernel = %s\nnGPinit = %d\nqBound = %s" % \
                    (nPar, sigma_d, whichOptim, tol_abs, kernelType, nGPinit,\
                     ', '.join(map(str, qBound))))

def nextGPsample(path2gpList):
    """
       Take the next sample of the parameters from their admissible space. 
       If the number of the available samples is less than a limit (=nGPinit), 
       take the initial samples randomly. Otherwise, use the BO-GP algorithm to draw the new sample. 
    """
    
    [xList,yList]=read_available_GPsamples(path2gpList,nPar)
    nData=len(yList)
    #xList=xList.reshape((nData,nPar))   #reshape as required by GPy and GPyOpt
    yList=yList.reshape((nData,1))       #reshape as required by GPy and GPyOpt
    
    if (nData<nGPinit):   #take initial random samples
       logger.info("take the sample randomly")
       tmp=[];
       for i in range(nPar):
          ##random initial sample
          minPar=domain[i]['domain'][0]
          maxPar=domain[i]['domain'][1]
          tmp.append(np.random.uniform(minPar,maxPar))
          ##some arbitrary value set by user
#           tmp.append(0.0)
       xNext=[];
       xNext.append(tmp)
    else:      #take GP samples  based on BO-GP algorithm
       #>>>> Define the GP model used in the BO
       ##kernel: (hyperparameters are optimizaed during the run)
       if kernelType=='Matern52':
          K=GPy.kern.Matern52(input_dim=nPar,lengthscale=1.0,variance=1.0)
       elif kernelType=='RBF':
          K=GPy.kern.RBF(input_dim=nPar,lengthscale=1.0,variance=1.0)

       gpModel = GPyOpt.models.gpmodel.GPModel(kernel=K,
                                        noise_var=sigma_d**2.,
                                        exact_feval=fevalFlag,
                                        optimizer='bfgs', #MLE optimization of the Kernel hyper parameters
                                        max_iters=200,
                                        optimize_restarts=5,
                                        verbose=False)

       #>>>> Set-up the BO (Bayesian Optimization) problem
       gprOpt=BayesianOptimization(f='',    #empty f
                               domain=domain,
                               maximize=maxFlag,   #default:False => minimization
                               model_type='GP',
                               model=gpModel,
                               acquisition_type='EI', # 'EI', 'MPI', 'UCB'
                               acquisition_jitter = 0.001,  #xi in acquisition func
                               acquisition_optimizer_type='lbfgs',  # to find max of acquistion func
                               exact_feval=fevalFlag,  #true: f_true is used with no noise
#                                inititial_design=init,   #random initial
                               X=xList,   #non-random initials
                               Y=ifac*yList,
                               normalize_Y=False,  #Normalize the outputs before performing any optimization.
                               verbosity=True)

       #Find the next x-sample
       xNext=gprOpt.suggest_next_locations(context=None, pending_X=None, ignored_X=None)
    logger.info('**** New GP sample is: %s' % ', '.join(map(str, xNext[0])))
    
    return np.array(xNext[0])

#///////////////////////////
def BO_update_convergence(xLast, yLast, path2gpList='./workDir/gpList.dat', \
                          path2figs='../figs/'):
    """
       1. Update gpList.dat by adding the last drawn sample and associated response to it
       2. Check if BO-GP is converged or not
       
       Note: exit(1) is taken as the signal of the convergence
    """
    [xList,yList]=read_available_GPsamples(path2gpList,nPar)
    nData=len(yList) # overwritten later
    #xList=xList.reshape((nData,nPar))   #reshape as required by GPy and GPyOpt
    yList=yList.reshape((nData,1))       #reshape as required by GPy and GPyOpt

    #Update gpList.dat
    update_GPsamples(path2gpList,xList,yList,xLast,yLast)
    #read the updated gpList.dat
    [xList_,yList_]=read_available_GPsamples(path2gpList,nPar)
    
    #Check convergence of BO
    #>>>>> plot convergence
    nData=len(yList_)
    if nData>1:
       [xDistList,yBestList]=\
           my_convergence_plot(xList_,yList_,whichOptim,path2figs,'bo_convergence')
      #>>>> Converged Optimal Value:
       # if (n>2): # check convergence
       #     err_d=xDistList[-1]
       #     if (err_d<tol_d):
       #         err_b=abs(yBestList[-1]-yBestList[-2])/abs(yBestList[-1])
       #         if (err_b<tol_b):
       #             xOpt=xList_[xList_.shape[0]-1]
       #             fxOpt=yList_[yList_.shape[0]-1]
       #             print(' ******* Converged Optimal Values (x,f(x))= (%g,%g)' % (xOpt,fxOpt))
       #             print('err_d, err_b=%f %f' %(err_d,err_b))
       #             #send convergence signal
       #             sys.exit(1)
    gpSurface_plot(xList_,yList_,nData)
    # check convergence: only for minimize !!!
    logger.debug("check convergence")
    if yList_[-1]<=tol_abs: # take into acc
        xOpt=xList_[-1]
        fxOpt=yList_[-1]
        logger.info(' ******* Converged Optimal Values q = %s, y = %f <= %f'\
                    % (', '.join(map(str, xOpt)), fxOpt, tol_abs))
        return 1
    else:
        logger.info("not converged yet, y = %f > %f" % (yList_[-1], tol_abs))
        return 0
    
#////////////////////////////////////
def gpSurface_plot(xList, yList, nData, path2figs="../figs/", Rmin=0, Rmax=0, \
                   bounds=BOUNDS):
    """ 
       Reconstruct the GPR and plot it in 2D (or 1D) planes of the parameters admissible space
    """
    #read the updated gpList.dat
    # [xList,yList]=read_available_GPsamples('./workDir/gpList.dat',nPar)
    # nData=len(yList)
    #xList=xList.reshape((nData,nPar))   #reshape as required by GPy and GPyOpt
    yList=yList.reshape((nData,1))       #reshape as required by GPy and GPyOpt

    #reshape the arrays according to GPy
    yGP=np.asarray(yList)
    if nPar==1:
        yGP=yGP[:,None]        #required by Gpy library
    xGP=np.asarray(xList)

    #>>>> Reconstruct and Plot GP surrogate in parameter space
    if nPar>1:
        #plot in 2D subspace of the parameters space
        plotOpts={'figDir':path2figs,
                  'figName':'gp2D_%02d'% (nData),
                  'kernelType':kernelType,   #required to construct the final GPR
                  'whichOptim':whichOptim}
        if not (Rmin==0 and Rmax==0):
            plotOpts["Rmin"]=Rmin
            plotOpts["Rmax"]=Rmax
        gpOpt2d_postProc(nPar,xGP,yGP,sigma_d,bounds,plotOpts)
    elif nPar==1:
        if Rmin==0 and Rmax==0:
            Rmin=0.05
            Rmax=10
        nTest=100   #no of test points, only for plot
        plotOpts={'figDir':path2figs,
                  'figName':'gp1D_%02d' % (nData),
                  'kernelType':kernelType,   #required to construct the final GPR
                  'whichOptim':whichOptim,
                  'arbitSample':'no',
                  'ylim':[Rmin,Rmax]} # NEED TO BE TUNED
        gpOpt1d_postProc(xGP,yGP,sigma_d,domain,nTest,plotOpts)
    else:
        logger.error("nPar should be >= 1")
        sys.exit(1)
