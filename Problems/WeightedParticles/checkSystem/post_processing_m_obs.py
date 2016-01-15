"""
Created on Wed Oct 28 10:05:37 2015

@author: pieter
"""
import time 
import tables
import scipy
import sys
import numpy as np
import matplotlib.pylab as plt

from scipy.linalg import norm
from scipy import sqrt

import matplotlib.pylab as plt


if __name__=="__main__":
    D = 1./2
    Dt = 1e-2    
    a = 1
    zeta = 0.
    alpha=1
    beta = -1
    lambd = scipy.array([a,D,alpha,beta,zeta])
   # param['eps'] = 1e-5
    #param['Dt'] = Dt
    #param['dt']=dt         


    rho_Dt_fp= np.loadtxt('data/25-11-Rho_pde_dx_e-4_eps-5.out') 
    Jv_pde= np.loadtxt('data/25-11-Jv_pde_dx_e-4_eps-5.out')
    
    sde_rho_sq = np.loadtxt('data/25-11-rho_sq_eps-5.out')
    E_rho = np.loadtxt('data/25-11-E_rho_eps-5.out')
    sq_E_rho  =np.loadtxt('data/25-11-sq_E_rho_eps-5.out')
    
    sde_Jv_sq = np.loadtxt('data/25-11-Jv_sq_eps-5.out')
    E_Jv = np.loadtxt('data/25-11-E_Jv_eps-5.out')   #lots of data, not needed when only interested in the variance
    sq_E_Jv =  np.loadtxt('data/25-11-sq_E_Jv_eps-5.out' )
    
    sde_Jv_sq_M10 = np.loadtxt('data/M10-Jv_sq_eps-5.out')[:-1]
    E_Jv_M10 = np.loadtxt('data/M10-E_Jv_eps-5.out')[:-1]   #lots of data, not needed when only interested in the variance
    sq_E_Jv_M10 =  np.loadtxt('data/M10-sq_E_Jv_eps-5.out' )[:-1]
            
        
    #SDE
    rho_Dt_sde = E_rho[-1]
    Jv_sde = E_Jv[-1]

    
    resize_factor = int (len(rho_Dt_fp)/len(rho_Dt_sde))
  #  resize_factor = 10
    print "Discretisation for solving sde is ",  resize_factor , " times coarser than the discretisation for solving the pde"
    
    rho_coarse = scipy.zeros(len(rho_Dt_sde))
    for i in range (0,len(rho_coarse)):
        bin_av = 0
        for j in range (0,  resize_factor ):
            bin_av = bin_av+ rho_Dt_fp[i*resize_factor+j]
        rho_coarse[i] = bin_av/resize_factor
        
    bins =  len(Jv_sde)
    Jv_coarse = scipy.zeros(bins)
    for i in range (0,len(Jv_coarse)):
        bin_av = 0
        for j in range (0,  resize_factor ):
            bin_av = bin_av+ Jv_pde[i*resize_factor+j]
        Jv_coarse[i] = bin_av/resize_factor
        
    pde_rho_sq = norm(rho_coarse)**2
    pde_Jv_sq = norm(Jv_coarse)**2

    
    
    log_flag = 'log'
  #  log_flag = 'linear'
    
   # Nlist = scipy.array([1000,2000,4000,8000,16000])
       
    Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000])
       
   # Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000, 128000])
#    plot_rms = plt.plot(Nlist, sde_rho_sq -pde_rho_sq, label = '$RMS^2$')
#    plt.xlabel('$N$', fontsize = 16)
#    plt.ylabel(r'MSE( \rho)$' , fontsize = 16)
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    plt.show()
#    plot_var = plt.plot(Nlist, sde_rho_sq - sq_E_rho, label = '$Variance$')
#    plt.ylabel(r'Var$(\rho)$' , fontsize = 16)
#    #plt.ylabel()
#    plt.xlabel('$N$', fontsize = 16)
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    plt.show()
#    plot_bias =  plt.plot(Nlist, abs(pde_rho_sq - sq_E_rho), label = '$Bias$')
#    plt.xlabel('$N$', fontsize = 16)
#    plt.ylabel(r'Bias $(\rho)$', fontsize = 16)
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    plt.show()
#    
    plot_rho = plt.plot( rho_coarse, label = '$\$')
    plot_rho =  plt.plot(rho_Dt_sde, label = '$\$')
    plt.ylabel(r' $(\rho)$', fontsize = 16)
 #   plt.xlabel('$x$', fontsize = 16)
    plt.show()
       
    plot_var = plt.plot(Nlist, (sde_Jv_sq - sq_E_Jv)/bins, label = '$\$')
    plot_var_M10 = plt.plot(Nlist, (sde_Jv_sq_M10 - sq_E_Jv_M10)/bins, label = '$\$')
    plt.ylabel('Var ($\mathbf{\hat{Jv}} $)' , fontsize = 16)
    plt.xlabel('$N$', fontsize = 16) 
    plt.xscale(log_flag)
    plt.yscale(log_flag)
    plt.show()
    
    bias =  scipy.zeros((len(Nlist), bins))
    biasnorm=  scipy.zeros(len(Nlist))
    bias_sq=  scipy.zeros(len(Nlist))
    
    for i in range (0,len(Nlist)):
        bias[i] = Jv_coarse - E_Jv[i]
     #   bias_M10[i] = Jv_coarse - E_Jv_M10[i]
        biasnorm[i] = norm(bias[i])
        bias_sq[i] = biasnorm[i]**2
        
    plot_biasnorm =  plt.plot(Nlist, biasnorm/bins, label = '$Biasnorm$')
    plot_biasnorm =  plt.plot(Nlist, biasnorm/bins, label = '$Biasnorm$')
    
    plt.xlabel('$N$', fontsize = 16)
    plt.ylabel(r'$\left(\mathbf{Bias}(\mathbf{\hat{Jv}} , \mathbf{\bar{Jv} }) \right)^2$' , fontsize = 16)
    
    
    
    
    
    plt.xscale(log_flag)
    plt.yscale(log_flag)
    
    
    
    plt.show()
    
    
    
    plot_Jv_mse = plt.plot(Nlist, (sde_Jv_sq - 2*np.dot(E_Jv,Jv_coarse)  + pde_Jv_sq)/bins, label = '$\$')
    plt.ylabel(r'MSE($\mathbf{\hat{Jv}} , \mathbf{\bar{Jv} }$)' , fontsize = 16)
    plt.xlabel('$N$', fontsize = 16)
    plt.xscale(log_flag)
    plt.yscale(log_flag)
    plt.show()
    
    plot_Jv_mse2 = plt.plot(Nlist,  (sde_Jv_sq - sq_E_Jv + bias_sq)/bins , label = '$\$')  #gives the smae result, as expected
    plt.ylabel(r'MSE($\mathbf{\hat{Jv}} , \mathbf{\bar{Jv} }$)' , fontsize = 16)
    plt.xlabel('$N$', fontsize = 16)
    plt.xscale(log_flag)
    plt.yscale(log_flag)
    plt.show()
    
    
    plot_Jv = plt.plot( Jv_coarse, label = '$\$')
    plot_Jv =  plt.plot(Jv_sde, label = '$\$')
    plt.ylabel(r' $Jv$', fontsize = 16)
   # plt.xlabel('$x$', fontsize = 16)
    plt.show()
#plt.xscale('log')
#plt.yscale('log')
#


#Jv_var = scipy.zeros(M)
#average= Jv_norm[-1]
#for m in range(1,M):
#    a= np.power((Jv_norm[m] - average),2) 
#    Jv_var[m]= np.sqrt(a)

#ax=plt.subplot(111)
#ax.set_xlim(1, 5000)
#ax.set_ylim(0, np.max(Jv_var))

#plot_var = plt.plot(Jv_var)
#plt.xlabel('$N$', fontsize = 16)
#plt.ylabel(r'$\frac{\mathbb{E}[ (\hat{\mathbf{\rho} } - \mathbf{\rho}_{FP} )^2 ] } {|| \mathbf{1}||^2}    $', fontsize = 28)
#plt.ylabel(r'$\mathbb{E}[ (\hat{\mathbf{\rho} } - \mathbf{\rho}_{FP} )^2 ]   $', fontsize = 18)
#plt.ylabel(r'$\mathbb{E}[ (\hat{\mathbf{J} } \cdot \mathbf{v}  - \mathbf{J_{FP}} \cdot \mathbf{v})^2 ]   $', fontsize = 18)
#plt.legend([plot5, plot6], loc='best')
#plt.legend(bbox_to_anchor=(1, 1), numpoints = 1 )
#plt.savefig("results/plots/ms_Jveps-5eps-6.pdf")
#plt.show()
                        
        

