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
from numpy import *

from scipy.linalg import norm
from scipy import sqrt

import matplotlib.pylab as plt


def get_fixpoint( mu, sigma, alpha):
    xi=0
    C=  sigma**2/(2*alpha)
    if (3*C>1):
        print "Make sure that 3D/(2alpha) <1 to get non-zero equilibrium solutions"
    else:
        xi = sqrt(1-3*C)*(1+ 6*mu/(sigma**2)*C**2*(1-2*C)/(1-3*C)) 
    return xi

if __name__=="__main__":
    sigma=1.
    mu = 0.1
    alpha= 1.5*sigma**2 + 1e-2

    xi_abs= get_fixpoint( mu, sigma, alpha)
    print 'Analytical solution = ' , xi_abs
       
   # Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000])
            
    xL = -1.7
    xR = 1.7
    dx = 1e-2
    grid = scipy.arange(xL+dx/2.,xR,dx)
    dt = 0.01
    
    Dt = 10000
   
    
    ######## Newton-bias
    
    xmean_list = zeros(len(NStates_1e6))
#    xmean_list_Dt1 = zeros(len(NStates_1e6))
#    xmean_list_Dt5 = zeros(len(NStates_1e6))
#    xmean_list_Dt10 = zeros(len(NStates_1e6))
#    xmean_list_Dt20 = zeros(len(NStates_1e6))
#    xmean_list_Dt50 = zeros(len(NStates_1e6))
#    
#    xmean_list_Dt1_1e5 = zeros(len(NStates_1e6))
#    xmean_list_Dt10_1e5 = zeros(len(NStates_1e6))
#    xmean_list_Dt1_1e4 = zeros(len(NStates_1e6))
#    xmean_list_Dt10_1e4 = zeros(len(NStates_1e6))
#    xmean_list_Dt10_1e7 = zeros(len(NStates_1e6))
#    xmean_list_Dt1_1e7 = zeros(len(NStates_1e6))
#        
        

    NStates_1e5_Dt10 =  np.loadtxt('Newton/critical_alpha_Newton_states_8it_N1e5_Dt_10')
    resnorm_1e5_Dt10 = np.loadtxt('Newton/critical_alpha_resnorm_Dt10_N1e5')
    
#    np.loadtxt('Newton/critical_alpha_Newton_states_8it_N1e6_Dt_10')
#    np.loadtxt('Newton/critical_alpha_resnorm_Dt10_N1e6')

   
    for k in range(len(NStates_1e5_Dt10)):
        xmean_list = np.dot(grid, NStates_1e5_Dt10[k])*dx     
#        xmean_list[k] = np.dot(grid, NStates_1e6[k])*dx  
#        xmean_list_Dt1[k] = np.dot(grid, NStates_1e6_Dt1[k])*dx  
#        xmean_list_Dt5[k] = np.dot(grid, NStates_1e6_Dt5[k])*dx 
#        xmean_list_Dt10[k] = np.dot(grid, NStates_1e6_Dt10[k])*dx  
#        xmean_list_Dt20[k] = np.dot(grid, NStates_1e6_Dt20[k])*dx     
#        xmean_list_Dt50[k] = np.dot(grid, NStates_1e6_Dt50[k])*dx     
#        
                 
        plt.plot(grid,NStates_1e5_Dt10[k], 'blue',  label=r'$\Delta T = 10$')
#        plt.plot(grid,NStates_1e6_Dt1[k], 'cyan', label=r'$\Delta T = 1$')  
#        plt.plot(grid,NStates_1e6_Dt5[k], 'red', label=r'$\Delta T = 5$')  
#        plt.plot(grid,NStates_1e6_Dt10[k], 'green', label=r'$\Delta T = 10$')
#        plt.plot(grid,NStates_1e6_Dt20[k], 'gray' , label=r'$\Delta T = 20$')
#        plt.plot(grid,NStates_1e6_Dt50[k], 'magenta' , label=r'$\Delta T = 50$')  
        plt.title(r' $k=%d$' %k)
        plot_anal =plt.axvline(- get_fixpoint(0.1, 1.0, alpha), color='red',  linestyle='--', linewidth=2, label = r'$\xi $' ) #${\bar{x}}_{analytical}$')
        plt.legend(prop={'size':10} )
#        if (k==1):
#            rho_ss_dir_DT50 = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_alpha5_mu0p05_Dt50')
#            plt.plot( grid, rho_ss_dir_DT50, 'black', label= r'Direct sim $\Delta T = 50$' )
#            plt.legend(prop={'size':10} )
#            plt.savefig('plots/Newton_states(DT)_it1_dirsim.pdf')
        plt.plot()
        plt.legend(prop={'size':10} )
      #  plt.plot(grid,rho_ss,'red' )
        plt.savefig('plots/critcial_Newton_states(DT)_it%d.pdf' %k)
        plt.show()
    
    plt.plot(abs(xmean_list + get_fixpoint(0.1, 1.0, alpha)), 'blue', label=r'$\Delta T = 10$')
#    plt.plot(abs(xmean_list_Dt1 + get_fixpoint(0.05, 1.0, 5.0)), 'cyan', label=r'$\Delta T = 1$')
#    plt.plot(abs(xmean_list_Dt5 + get_fixpoint(0.05, 1.0, 5.0)), 'red', label=r'$\Delta T = 5$')
#    plt.plot(abs(xmean_list_Dt10 + get_fixpoint(0.05, 1.0, 5.0)), 'green', label=r'$\Delta T = 10$')
#    plt.plot(abs(xmean_list_Dt20 + get_fixpoint(0.05, 1.0, 5.0)), 'gray', label=r'$\Delta T = 20$')
#    plt.plot(abs(xmean_list_Dt50 + get_fixpoint(0.05, 1.0, 5.0)), 'magenta', label=r'$\Delta T = 50$')   
 #   plt.legend(bbox_to_anchor=(1, 0.64), numpoints = 1,  prop={'size':10} )
    plt.ylabel(r'$ |\bar{x} - \xi |$', size=15)
    plt.xlabel('k')
    plt.title(r'$N=10^6$')
   # plt.ylim(0,0.6)
    plt.savefig('plots/critical_Newton_states_abs_bias_N1e6.pdf')
  #  plt.yscale('log')
    
    plt.show()
#    DtList = [0.1, 1, 2, 10, 20, 50, 100, 200, 1000]
#    Biaslist = [ np.dot(grid, evolved_state0)*dx , np.dot(grid, evolved_state1)*dx , np.dot(grid, evolved_state2)*dx , np.dot(grid, evolved_state10)*dx ,
#                np.dot(grid, evolved_state20)*dx , np.dot(grid, evolved_state50)*dx , np.dot(grid, evolved_state100)*dx ,np.dot(grid, evolved_state200)*dx ,np.dot(grid, evolved_state1000)*dx   ] 
#    Biaslist = Biaslist + get_fixpoint(0.05, 1.0, 5.0)
#
#    plt.xscale('log')   
#    plt.yscale('log')
#    plt.plot(DtList, Biaslist, label='direct simulation', linewidth=1.5)
#    Dtlist_Newton_1 = np.array([0.1, 1, 5, 10, 20, 50])
#    #Newtonbias_list_8 = [xmean_list[-1] , xmean_list_Dt1[-1] , xmean_list_Dt5[-1] , xmean_list_Dt10[-1] , xmean_list_Dt20[-1] , xmean_list_Dt50[-1] ] 
#    Newtonbias_list = np.array([xmean_list , xmean_list_Dt1 , xmean_list_Dt5 , xmean_list_Dt10 , xmean_list_Dt20 , xmean_list_Dt50 ] )
#    Newtonbias_list_T = np.transpose(Newtonbias_list )
#    
#    plt.xlabel(r'$\Delta T$')
#    plt.ylabel(r'$ \bar{x} - \xi $', size=16)
#    for k in range(1,9):
#            plt.plot(Dtlist_Newton_1*k, Newtonbias_list_T[k] + get_fixpoint(0.05, 1.0, 5.0) , label='k=%d' %k )
#           # plt.plot(Dtlist_Newton_1, xmean_list_Dt50[k] + get_fixpoint(0.05, 1.0, 5.0), 'red')
#    plt.legend(bbox_to_anchor=(0.4,0.8), numpoints = 1,  prop={'size':10} )
#    plt.savefig('plots/bias_compare_Newton_direct.pdf')
#    plt.show()
#    
       ######## Newton-bias
    plt.plot(abs(xmean_list + get_fixpoint(0.1, 1.0, alpha)), 'red', linestyle= '-' , label=r'$N=10^7$')
    #plt.text(3, 0.09, r'$\Delta T=10$')
    #plt.text(3, 0.52, r'$\Delta T=1$')
    plt.legend(bbox_to_anchor=(1, 0.64), numpoints = 1,  prop={'size':10} )
    plt.ylabel(r'$ |\bar{x} - \xi| $', size=16)
    plt.xlabel('k')
   # plt.ylim(-0.05,0.6)
    plt.savefig('plots/Critical_Newton_states_abs_bias_N1e5.pdf')
    
    plt.show()

    
    ###########    RESIDUAL         ######################################################

    resnorm_N1e5_Dt10 = np.loadtxt('Newton/critical_alpha_resnorm_Dt10_N1e5')
   
    plt.yscale('log')
    plt.plot(     resnorm_N1e5_Dt10  ,  "ko-" , label = r'$N=10^5$',markersize=4)    
#    plt.plot(    resnorm_DT10_N1e3 ,  "yo-" , label = r'$N=10^3$',markersize=4)    
#    plt.plot(    resnorm_DT10_N1e4 ,  "bo-" , label = r'$N=10^4$',markersize=4)       
#    plt.plot(    resnorm_DT10_N1e5 ,  "ro-" , label = r'$N=10^5$',markersize=4)    
#    plt.plot(    resnorm_DT10_N1e6 ,  "co-" , label = r'$N=10^6$',markersize=4)    
#    plt.plot(    resnorm_DT10_N1e7 ,  "mo-" , label = r'$N=10^7$',markersize=4) 
    plt.legend(  prop={'size':9} )
    plt.xlabel(r'$k$', fontsize = 14)
    plt.ylabel(r'$||  \rho - \Phi^N_T(\rho)||$'  , fontsize = 14)
    plt.savefig("plots/critical_res_norm_gmres_e-5_alpha5_mu0p05.pdf")
    plt.show()
    

    
