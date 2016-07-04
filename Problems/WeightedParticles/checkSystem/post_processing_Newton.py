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

def get_density( mu, sigma, alpha, grid):
    xi=0
    dx= grid[-1]-grid[-2]
    C=  sigma**2/(2*alpha)
    if (3*C>1):
        print "Make sure that 3D/(2alpha) <1 to get non-zero equilibrium solutions"
    else:
        xi = sqrt(1-3*C)*(1+ 6*mu/(sigma**2)*C**2*(1-2*C)/(1-3*C)) 
          
    norm_c = sum( np.exp(-alpha*(grid-xi)**2/sigma**2 + (-grid**4 + grid**2)*2*mu/(sigma**2) ))*dx
    rho_ss1 = np.exp(-alpha*(grid-xi)**2/sigma**2 + (-grid**4 + grid**2)*2*mu/(sigma**2))/norm_c
    xi=-xi
    rho_ss2 =  np.exp(-alpha*(grid-xi)**2/sigma**2 + (-grid**4 + grid**2)*2*mu/(sigma**2))/norm_c
    return rho_ss1, rho_ss2


if __name__=="__main__":
    sigma=1.
    mu = 0.05
    alpha= 5.

    xi_abs= get_fixpoint( mu, sigma, alpha)
    print 'Analytical solution = ' , xi_abs
       
   # Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000])
       
    xL = -1.7
    xR = 1.7
    dx = 1e-2
    grid = scipy.arange(xL+dx/2.,xR,dx)
    dt = 0.01
    
    Dt = 10000
    timesteps = scipy.arange(dt,Dt+dt,dt)
    xi1 = full(Dt/dt, xi_abs)
    xi2 = full( Dt/dt, -xi_abs)
    dt = 1e-2 # ti
    
    alphas = scipy.arange(0, 6, 0.11)
    fixpoints = scipy.zeros_like(alphas)
    for i in range (len(fixpoints)):
        fixpoints[i] = get_fixpoint(mu, sigma, alphas[i])
    
    plot_bif = plt.plot(alphas, fixpoints)
    plot_bif = plt.plot(alphas, -fixpoints)
    plt.ylabel(r'$\xi$', fontsize = 18)
    plt.xlabel(r'$\alpha$', fontsize = 16)
  #  plt.savefig('plots/bif_alpha.pdf')
    plt.show()

    rho_ss_alpha1= np.loadtxt('pfix_alpha1.txt')  
    xmean_alpha1 = np.dot(grid,rho_ss_alpha1)*dx         
    rho_ss_alpha5 =np.loadtxt('pfix_alpha5.txt')
    xmean_alpha5 = np.dot(grid,rho_ss_alpha5)*dx
    plt.axvline(xmean_alpha1, linewidth=2)
    plt.plot()
    plot_anal =plt.axvline( get_fixpoint( mu, sigma, 1), color='red',  linestyle='--', linewidth=2, label =  r'${\xi}$')
    rho_ss_anal_1 = get_density(mu, sigma, 1, grid)[1]
    plot_anal_density = plt.plot(grid, rho_ss_anal_1 , color='orange' , label = r'$\rho_{\xi}(x)$')
    
    

    plot_stable = plt.plot(grid, rho_ss_alpha1, label = r'${\rho^*}_{\alpha=1}$')
    plt.axvline(xmean_alpha5, color='green', linewidth=2)
    plt.axvline( -get_fixpoint( mu, sigma, 5), color='red',  linestyle='--', linewidth=2)
    rho_ss_anal_5 = get_density(mu, sigma, 5, grid)[1]
    plot_anal_density = plt.plot(grid,  rho_ss_anal_5 , color='orange' , label = r'$\rho_{\xi}(x)$')
    plt.xlabel('$x$', fontsize = 16)
    plt.ylabel(r'$\rho*$', fontsize = 18)
    delta_x_alpha1 = np.abs(xmean_alpha1- get_fixpoint( mu, sigma, 1))
    delta_x_alpha5  =  np.abs(xmean_alpha5 + xi_abs)
    #plt.title( r'$\Delta \bar{x} (\alpha=5) =%.3f}$    ' %delta_x_alpha5 + r'$\Delta \bar{x} (\alpha=1) =%.3f}$' %delta_x_alpha1 , fontsize=11)
    plot_meta= plt.plot(grid, rho_ss_alpha5, label = r'${\rho^*}_{\alpha=5}$')
    plt.legend([plot_anal, plot_meta, plot_stable], loc='best')

   # plt.legend(bbox_to_anchor=(1, 1), numpoints = 1 )len(NStates_1e6
    plt.legend(loc=0, prop={'size':11}) #fontsize = 'x-small')
    plt.gca().set_ylim([0, 1.7])   
#    plt.savefig('plots/steady_states.pdf')
    plt.savefig('plots/steady_states_upd2.pdf')
    plt.show()    
    
    print 'Analytical solution: systemic risk = ', xi_abs
    print 'Newton-Krylov solution: systemic risk = ', xmean_alpha5
  
   # plot_sysrisk =  plt.plot( timesteps, xmean  , label = '$\$')
    plot_sysrisk = plt.plot( timesteps, xi1, linestyle='--', linewidth = 2 , color = 'red', label = '$\$')
    plot_sysrisk = plt.plot( timesteps, xi2, linestyle='--', linewidth = 2,  color = 'red' , label = '$\$')
    #plot_rho =  plt.plot(rho_Dt_sde, label = '$\$')
    plt.ylabel(r' $\bar{x}$', fontsize = 18)
    plt.xlabel('$t$', fontsize = 16)
    plt.title(r'$N=100$, $\mu=%.1f$, $\sigma=%d$, $\alpha=%d$' %(mu, sigma, alpha) )
 #  plt.savefig('plots/xmean_alpha%d.pdf'%alpha, bbox_inches='tight')
    plt.show()
    
    
    
    #    DIRECT SIMULATION  #######################################
    
    NStates_1e6 = np.loadtxt('Newton/Newton_states_8it_N1e6_Dt1e-1')
    NStates_1e6_Dt1 = np.loadtxt('Newton/Newton_states_8it_N1e6_Dt_1')
    NStates_1e6_Dt5 = np.loadtxt('Newton/Newton_states_8it_N1e6_Dt_5')
    NStates_1e6_Dt10 = np.loadtxt('Newton/Newton_states_8it_N1e6_Dt_10')
    NStates_1e6_Dt20 = np.loadtxt('Newton/Newton_states_8it_N1e6_Dt_20')
    NStates_1e6_Dt50 = np.loadtxt('Newton/Newton_states_8it_N1e6_Dt_50')
    evolved_state_obs = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_Dt2_alpha5_diff_force')    
   # evolved_state = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_Dt2_alpha10_mu2')   
    evolved_state0 = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_alpha5_mu0p05_Dt0')    #0.1
    evolved_state1 = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_alpha5_mu0p05_Dt1')  
    evolved_state2 = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_Dt2_alpha5_mu0p05')
    evolved_state10 = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_alpha5_mu0p05_Dt10')  
    evolved_state20 = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_alpha5_mu0p05_Dt20')  
    evolved_state50 = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_alpha5_mu0p05_Dt50')  
    evolved_state100 = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_alpha5_mu0p05_Dt100')  
    evolved_state200 = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_alpha5_mu0p05_Dt200')  
    evolved_state1000 = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_alpha5_mu0p05_Dt1000')  
    #evolved_state = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_Dt2_alpha5_mu0p05')  


    plt.plot(grid,evolved_state0, 'gray', label=r'$\Delta T = 0.1$')
    plt.plot(grid,evolved_state1, 'cyan' , label=r'$\Delta T = 1$')
    plt.plot(grid,evolved_state2, 'red', linestyle='--',  label=r'$\Delta T = 2$')
    plt.plot(grid,evolved_state10, 'black', label=r'$\Delta T = 10$')
    plt.plot(grid,evolved_state20, 'orange',  label=r'$\Delta T = 20$')
    plt.plot(grid,evolved_state50, 'blue', linestyle='--',  label=r'$\Delta T = 50$')
    plt.plot(grid,evolved_state100, 'green',  label=r'$\Delta T = 100$')
    plt.plot(grid,evolved_state200, 'blue',  label=r'$\Delta T =200$')
    plt.plot(grid,evolved_state1000, 'yellow', linestyle='--',  label=r'$\Delta T =1000$')
    plt.plot(grid, NStates_1e6[-1], 'magenta', linewidth=4,  label=r'Newton')
    plt.plot(grid, NStates_1e6_Dt10[-1], 'magenta', linewidth=4,  label=r'Newton')
   #     plt.axvline(xmean, linewidth=2, color='gray')
    plt.axvline(np.dot(grid, evolved_state0)*dx   , linewidth=1, color='gray')  
    plt.axvline(np.dot(grid, evolved_state1)*dx   , linewidth=1, color='cyan')
  #  plt.axvline(np.dot(grid, evolved_state2)*dx   , linewidth=1, color='red')
    plt.axvline(np.dot(grid, evolved_state10)*dx   , linewidth=1, color='black')
    plt.axvline(np.dot(grid, evolved_state100)*dx   , linewidth=2, color='green')    
    plt.axvline(np.dot(grid, evolved_state200)*dx   , linewidth=2, color='blue')    
    plot_anal =plt.axvline(- get_fixpoint(0.05, 1.0, 5.0), color='red',  linestyle='--', linewidth=2, label =  r'${\bar{x}}_{analytical}$')

    
    plt.legend()
 #   plt.savefig('plots/steady_states_f(Dt)_obs.pdf')
    plt.show()
    
    
    plt.plot(grid,evolved_state1, 'cyan' , label=r'$\Delta T = 1$')
    plt.plot(grid,evolved_state10, 'black', label=r'$\Delta T = 10$')
    plt.plot(grid,evolved_state20, 'orange',  label=r'$\Delta T = 20$')
    plt.plot(grid,evolved_state50, 'blue', linestyle='--',  label=r'$\Delta T = 50$')
    plt.plot(grid,evolved_state100, 'green',  label=r'$\Delta T =100$')
    plt.plot(grid,evolved_state1000, 'orange', linestyle='--',  label=r'$\Delta T =1000$')
   #     plt.axvline(xmean, linewidth=2, color='gray')
    plt.axvline(np.dot(grid, evolved_state1)*dx   , linewidth=1, color='cyan')
    plt.axvline(np.dot(grid, evolved_state10)*dx   , linewidth=1, color='black')
    plt.axvline(np.dot(grid, evolved_state100)*dx   , linewidth=2, color='green')    
    plt.axvline(np.dot(grid, evolved_state1000)*dx   , linewidth=2, color='blue')    
    plot_anal =plt.axvline(- get_fixpoint(0.05, 1.0, 5.0), color='red',  linestyle='--', linewidth=2, label =  r'$\xi$')
    
    
    plt.legend()
    plt.savefig('plots/steady_states_f(Dt).pdf')
    plt.show()
    
    
    DtList = [0.1, 1, 2, 10, 20, 50, 100, 200, 1000]
    Biaslist = [ np.dot(grid, evolved_state0)*dx , np.dot(grid, evolved_state1)*dx , np.dot(grid, evolved_state2)*dx , np.dot(grid, evolved_state10)*dx ,
                np.dot(grid, evolved_state20)*dx , np.dot(grid, evolved_state50)*dx , np.dot(grid, evolved_state100)*dx ,np.dot(grid, evolved_state200)*dx ,np.dot(grid, evolved_state1000)*dx   ] 
    Biaslist = Biaslist + get_fixpoint(0.05, 1.0, 5.0)
    
    plt.xscale('log')   
    plt.yscale('log')
    plt.plot(DtList, Biaslist)
    plt.xlabel(r'$\Delta T$')
    plt.ylabel(r'$ \bar{x} - \xi $', size=16)
    plt.savefig('plots/bias_steady_states_f(Dt).pdf')
    plt.show()
    
    #Time evolution of mean during one simulation     
    x_mean200 = np.loadtxt('Newton/x_mean_N1e6_mu0p05_lambda5_Dt200')
    x_mean200_N1e5  = np.loadtxt('Newton/x_mean_N1e5_mu0p05_alpha5_Dt200')
    plt.ylabel(r'$| \bar{x}- \xi |$', size=16)
    plt.xlabel(r'$\Delta t$')
    plt.plot(abs(x_mean200+ get_fixpoint(0.05, 1.0, 5.0)), label=r'$N=10^6$')
    plt.plot(abs(x_mean200_N1e5 + get_fixpoint(0.05, 1.0, 5.0)), label=r'$N=10^5$')
    
#    plt.xscale('log')
    plt.legend()
   # plt.yscale('log')
    plt.savefig('plots/bias_steady_states_f(dt)_N1e5_6.pdf')
    plt.show()
    
    x_mean20 = np.loadtxt('Newton/x_mean_N1e7_mu0p05_alpha5_Dt20')
    plt.ylabel(r'$ \bar{x}- \xi $', size=16)
    plt.xlabel(r'$\Delta t$')
    plt.plot(x_mean20+ get_fixpoint(0.05, 1.0, 5.0))
    plt.savefig('plots/bias_steady_states_f(dt)_N1e7.pdf')
    plt.show()
    
   # 
    
    ######## Newton-bias
    
    xmean_list = zeros(len(NStates_1e6))
    xmean_list_Dt1 = zeros(len(NStates_1e6))
    xmean_list_Dt5 = zeros(len(NStates_1e6))
    xmean_list_Dt10 = zeros(len(NStates_1e6))
    xmean_list_Dt20 = zeros(len(NStates_1e6))
    xmean_list_Dt50 = zeros(len(NStates_1e6))
    
    xmean_list_Dt1_1e5 = zeros(len(NStates_1e6))
    xmean_list_Dt10_1e5 = zeros(len(NStates_1e6))
    xmean_list_Dt1_1e4 = zeros(len(NStates_1e6))
    xmean_list_Dt10_1e4 = zeros(len(NStates_1e6))
    xmean_list_Dt10_1e7 = zeros(len(NStates_1e6))
    xmean_list_Dt1_1e7 = zeros(len(NStates_1e6))
        
    NStates_1e7_Dt1 = np.loadtxt('Newton/Newton_states_8it_N1e7_Dt_1')
    NStates_1e7_Dt10 = np.loadtxt('Newton/Newton_states_8it_N1e7_Dt_10')
    NStates_1e5_Dt1 = np.loadtxt('Newton/Newton_states_8it_N1e5_Dt_1')
    NStates_1e5_Dt10 = np.loadtxt('Newton/Newton_states_8it_N1e5_Dt_10')
    NStates_1e4_Dt1 = np.loadtxt('Newton/Newton_states_8it_N1e4_Dt_1')
    NStates_1e4_Dt10 = np.loadtxt('Newton/Newton_states_8it_N1e4_Dt_10')
    
    for k in range(len(NStates_1e6)):
        xmean = np.dot(grid, NStates_1e6[k])*dx     
        xmean_list[k] = np.dot(grid, NStates_1e6[k])*dx  
        xmean_list_Dt1[k] = np.dot(grid, NStates_1e6_Dt1[k])*dx  
        xmean_list_Dt5[k] = np.dot(grid, NStates_1e6_Dt5[k])*dx 
        xmean_list_Dt10[k] = np.dot(grid, NStates_1e6_Dt10[k])*dx  
        xmean_list_Dt20[k] = np.dot(grid, NStates_1e6_Dt20[k])*dx     
        xmean_list_Dt50[k] = np.dot(grid, NStates_1e6_Dt50[k])*dx     
        
        xmean_list_Dt1_1e7[k] = np.dot(grid, NStates_1e7_Dt1[k])*dx     
        xmean_list_Dt10_1e7[k] = np.dot(grid, NStates_1e7_Dt10[k])*dx     
        xmean_list_Dt1_1e5[k] = np.dot(grid, NStates_1e5_Dt1[k])*dx  
        xmean_list_Dt10_1e5[k] = np.dot(grid, NStates_1e5_Dt10[k])*dx  
        xmean_list_Dt1_1e4[k] = np.dot(grid, NStates_1e4_Dt1[k])*dx  
        xmean_list_Dt10_1e4[k] = np.dot(grid, NStates_1e4_Dt10[k])*dx 
        #xmean_evol = np.dot(grid, evolved_state0)*dx     
        plt.plot(grid,NStates_1e6[k], 'blue',  label=r'$\Delta T = 0.1$')
        plt.plot(grid,NStates_1e6_Dt1[k], 'cyan', label=r'$\Delta T = 1$')  
        plt.plot(grid,NStates_1e6_Dt5[k], 'red', label=r'$\Delta T = 5$')  
        plt.plot(grid,NStates_1e6_Dt10[k], 'green', label=r'$\Delta T = 10$')
        plt.plot(grid,NStates_1e6_Dt20[k], 'gray' , label=r'$\Delta T = 20$')
        plt.plot(grid,NStates_1e6_Dt50[k], 'magenta' , label=r'$\Delta T = 50$')  
        plt.title(r' $k=%d$' %k)
        plot_anal =plt.axvline(- get_fixpoint(0.05, 1.0, 5.0), color='black',  linestyle='--', linewidth=2, label = r'$\xi $' ) #${\bar{x}}_{analytical}$')
        plt.legend(prop={'size':10} )
#        if (k==1):
#            rho_ss_dir_DT50 = np.loadtxt('Newton/Newton_states_rho_ss_N1e6_alpha5_mu0p05_Dt50')
#            plt.plot( grid, rho_ss_dir_DT50, 'black', label= r'Direct sim $\Delta T = 50$' )
#            plt.legend(prop={'size':10} )
#            plt.savefig('plots/Newton_states(DT)_it1_dirsim.pdf')
        plt.plot()
        plt.legend(prop={'size':10} )
      #  plt.plot(grid,rho_ss,'red' )
        
        rho_ss_anal_5 = get_density(0.05, 1.0 , 5, grid)[1]
        plot_anal_density = plt.plot(grid,  rho_ss_anal_5 , color='black' , linestyle = '--' , label = r'$\rho_{\xi}(x)$')
        plt.savefig('plots/Newton_states(DT)_it%d_with_analytic_rho.pdf' %k)
        plt.show()
    
    plt.plot(abs(xmean_list + get_fixpoint(0.05, 1.0, 5.0)), 'blue', label=r'$\Delta T = 0.1$')
    plt.plot(abs(xmean_list_Dt1 + get_fixpoint(0.05, 1.0, 5.0)), 'cyan', label=r'$\Delta T = 1$')
    plt.plot(abs(xmean_list_Dt5 + get_fixpoint(0.05, 1.0, 5.0)), 'red', label=r'$\Delta T = 5$')
    plt.plot(abs(xmean_list_Dt10 + get_fixpoint(0.05, 1.0, 5.0)), 'green', label=r'$\Delta T = 10$')
    plt.plot(abs(xmean_list_Dt20 + get_fixpoint(0.05, 1.0, 5.0)), 'gray', label=r'$\Delta T = 20$')
    plt.plot(abs(xmean_list_Dt50 + get_fixpoint(0.05, 1.0, 5.0)), 'magenta', label=r'$\Delta T = 50$')   
    plt.legend(bbox_to_anchor=(1, 0.64), numpoints = 1,  prop={'size':10} )
    plt.ylabel(r'$ |\bar{x} - \xi |$', size=15)
    plt.xlabel('k')
    plt.title(r'$N=10^6$')
   # plt.ylim(0,0.6)
    plt.savefig('plots/Newton_states_abs_bias_N1e6.pdf')
  #  plt.yscale('log')
    
    plt.show()
    DtList = [0.1, 1, 2, 10, 20, 50, 100, 200, 1000]
    Biaslist = [ np.dot(grid, evolved_state0)*dx , np.dot(grid, evolved_state1)*dx , np.dot(grid, evolved_state2)*dx , np.dot(grid, evolved_state10)*dx ,
                np.dot(grid, evolved_state20)*dx , np.dot(grid, evolved_state50)*dx , np.dot(grid, evolved_state100)*dx ,np.dot(grid, evolved_state200)*dx ,np.dot(grid, evolved_state1000)*dx   ] 
    Biaslist = Biaslist + get_fixpoint(0.05, 1.0, 5.0)

    plt.xscale('log')   
    plt.yscale('log')
    plt.plot(DtList, Biaslist, label='direct simulation', linewidth=1.5)
    Dtlist_Newton_1 = np.array([0.1, 1, 5, 10, 20, 50])
    #Newtonbias_list_8 = [xmean_list[-1] , xmean_list_Dt1[-1] , xmean_list_Dt5[-1] , xmean_list_Dt10[-1] , xmean_list_Dt20[-1] , xmean_list_Dt50[-1] ] 
    Newtonbias_list = np.array([xmean_list , xmean_list_Dt1 , xmean_list_Dt5 , xmean_list_Dt10 , xmean_list_Dt20 , xmean_list_Dt50 ] )
    Newtonbias_list_T = np.transpose(Newtonbias_list )
    
    plt.xlabel(r'$\Delta T$')
    plt.ylabel(r'$ \bar{x} - \xi $', size=16)
    for k in range(1,9):
            plt.plot(Dtlist_Newton_1*k, Newtonbias_list_T[k] + get_fixpoint(0.05, 1.0, 5.0) , label='k=%d' %k )
           # plt.plot(Dtlist_Newton_1, xmean_list_Dt50[k] + get_fixpoint(0.05, 1.0, 5.0), 'red')
    plt.legend(bbox_to_anchor=(0.4,0.8), numpoints = 1,  prop={'size':10} )
    plt.savefig('plots/bias_compare_Newton_direct.pdf')
    plt.show()
    
       ######## Newton-bias
    plt.plot(abs(xmean_list_Dt10_1e7 + get_fixpoint(0.05, 1.0, 5.0)), 'red', linestyle= '-' , label=r'$N=10^7$')
    plt.plot(abs(xmean_list_Dt1_1e7 + get_fixpoint(0.05, 1.0, 5.0)), 'red' , linestyle= '-' )#label=r'$\Delta T = 1, N=10^7$')
    plt.plot(abs(xmean_list_Dt1 + get_fixpoint(0.05, 1.0, 5.0)), 'green', linestyle= '--',  label=r'$N=10^6$')
    plt.plot(abs(xmean_list_Dt10 + get_fixpoint(0.05, 1.0, 5.0)), 'green', linestyle= '--'  )
    plt.plot(abs(xmean_list_Dt1_1e5 + get_fixpoint(0.05, 1.0, 5.0)), 'blue' , linestyle= ':', label=r'$N=10^5$')
    plt.plot(abs(xmean_list_Dt10_1e5 + get_fixpoint(0.05, 1.0, 5.0)), 'blue' , linestyle= ':')
    plt.plot(abs(xmean_list_Dt1_1e4 + get_fixpoint(0.05, 1.0, 5.0)), 'cyan',  label=r'$N=10^4$')
    plt.plot(abs(xmean_list_Dt10_1e4 + get_fixpoint(0.05, 1.0, 5.0)), 'cyan') # linestyle= '-.' )   
 #   plt.title(r'$\Delta T= 1$, $\Delta T=10$')
    plt.text(3, 0.09, r'$\Delta T=10$')
    plt.text(3, 0.52, r'$\Delta T=1$')
    plt.legend(bbox_to_anchor=(1, 0.64), numpoints = 1,  prop={'size':10} )
    plt.ylabel(r'$ |\bar{x} - \xi| $', size=16)
    plt.xlabel('k')
   # plt.ylim(-0.05,0.6)
    plt.savefig('plots/Newton_states_abs_bias_N1e4_5_6_7.pdf')
    
    plt.show()
    
    ###Log_plots
    
    plt.plot(abs(xmean_list + get_fixpoint(0.05, 1.0, 5.0)), 'blue', label=r'$\Delta T = 0.1$')
    plt.plot(abs(xmean_list_Dt1 + get_fixpoint(0.05, 1.0, 5.0)), 'cyan', label=r'$\Delta T = 1$')
    plt.plot(abs(xmean_list_Dt5 + get_fixpoint(0.05, 1.0, 5.0)), 'red', label=r'$\Delta T = 5$')
    plt.plot(abs(xmean_list_Dt10 + get_fixpoint(0.05, 1.0, 5.0)), 'green', label=r'$\Delta T = 10$')
    plt.plot(abs(xmean_list_Dt20 + get_fixpoint(0.05, 1.0, 5.0)), 'gray', label=r'$\Delta T = 20$')
    plt.plot(abs(xmean_list_Dt50 + get_fixpoint(0.05, 1.0, 5.0)), 'magenta', label=r'$\Delta T = 50$')   
    plt.legend(bbox_to_anchor=(0.8, 0.84), numpoints = 1,  prop={'size':10} )
    plt.ylabel(r'$ |\bar{x} - \xi |$', size=15)
    plt.xlabel('k')
    plt.title(r'$N=10^6$')
   # plt.ylim(0,0.6)
    plt.yscale('log')
    
    plt.savefig('plots/Newton_states_abs_bias_N1e6_log.pdf')

    plt.show()
    
       ######## Newton-bias
    
    plt.plot(abs(xmean_list_Dt10_1e7 + get_fixpoint(0.05, 1.0, 5.0)), 'red', linestyle= '-' , label=r'$N=10^7$')
    plt.plot(abs(xmean_list_Dt1_1e7 + get_fixpoint(0.05, 1.0, 5.0)), 'red' , linestyle= '-' )#label=r'$\Delta T = 1, N=10^7$')
    plt.plot(abs(xmean_list_Dt1 + get_fixpoint(0.05, 1.0, 5.0)), 'green', linestyle= '--',  label=r'$N=10^6$')
    plt.plot(abs(xmean_list_Dt10 + get_fixpoint(0.05, 1.0, 5.0)), 'green', linestyle= '--'  )
    plt.plot(abs(xmean_list_Dt1_1e5 + get_fixpoint(0.05, 1.0, 5.0)), 'blue' , linestyle= ':', label=r'$N=10^5$')
    plt.plot(abs(xmean_list_Dt10_1e5 + get_fixpoint(0.05, 1.0, 5.0)), 'blue' , linestyle= ':')
    plt.plot(abs(xmean_list_Dt1_1e4 + get_fixpoint(0.05, 1.0, 5.0)), 'cyan',  label=r'$N=10^4$')
    plt.plot(abs(xmean_list_Dt10_1e4 + get_fixpoint(0.05, 1.0, 5.0)), 'cyan') # linestyle= '-.' )   

 
    plt.legend(bbox_to_anchor=(1, 0.84), numpoints = 1,  prop={'size':10} )
    plt.ylabel(r'$ |\bar{x} - \xi| $', size=16)
    plt.xlabel('k')
   # plt.ylim(-0.05,0.6)
    plt.yscale('log')
    plt.savefig('plots/Newton_states_abs_bias_N1e4_5_6_7_log.pdf')
    

    
    plt.show()
    ###########    RESIDUAL         ######################################################

    resnorm_DT10_N1e2 = np.loadtxt('Newton/new_method_resnorm_Dt10_N1e2')
    resnorm_DT10_N1e3 = np.loadtxt('Newton/new_method_resnorm_Dt10_N1e3')
    resnorm_DT10_N1e4 = np.loadtxt('Newton/new_method_resnorm_Dt10_N1e4')
    resnorm_DT10_N1e5 = np.loadtxt('Newton/new_method_resnorm_Dt10_N1e5')
    resnorm_DT10_N1e6 = np.loadtxt('Newton/new_method_resnorm_Dt10_N1e6')
    resnorm_DT10_N1e7 = np.loadtxt('Newton/new_method_resnorm_Dt10_N1e7')  
   
    plt.yscale('log')
    plt.plot(    resnorm_DT10_N1e2 ,  "ko-" , label = r'$N=10^3$',markersize=4)    
    plt.plot(    resnorm_DT10_N1e3 ,  "yo-" , label = r'$N=10^3$',markersize=4)    
    plt.plot(    resnorm_DT10_N1e4 ,  "bo-" , label = r'$N=10^4$',markersize=4)       
    plt.plot(    resnorm_DT10_N1e5 ,  "ro-" , label = r'$N=10^5$',markersize=4)    
    plt.plot(    resnorm_DT10_N1e6 ,  "co-" , label = r'$N=10^6$',markersize=4)    
    plt.plot(    resnorm_DT10_N1e7 ,  "mo-" , label = r'$N=10^7$',markersize=4) 
    plt.legend(  prop={'size':9} )
    plt.xlabel(r'$k$', fontsize = 14)
    plt.ylabel(r'$||  \rho - \Phi^N_T(\rho)||$'  , fontsize = 14)
    plt.savefig("plots/res_norm_gmres_e-5_alpha5_mu0p05.pdf")
    plt.show()
    
    ##CONVERGENCE OF BIAS
    res_norm_mean_N1e2 = sum(resnorm_DT10_N1e2[5:])/len(resnorm_DT10_N1e3[5:]) 
    res_norm_mean_N1e3 = sum(resnorm_DT10_N1e3[5:])/len(resnorm_DT10_N1e3[5:]) 
    res_norm_mean_N1e4= sum(resnorm_DT10_N1e4[5:])/len(resnorm_DT10_N1e4[5:]) 
    res_norm_mean_N1e5 = sum(resnorm_DT10_N1e5[5:])/len(resnorm_DT10_N1e5[5:]) 
    res_norm_mean_N1e6 = sum(resnorm_DT10_N1e6[5:])/len(resnorm_DT10_N1e6[5:])
    res_norm_mean_N1e7 = sum(resnorm_DT10_N1e7[5:])/len(resnorm_DT10_N1e7[5:])
    
    res_norm_N = np.array([     res_norm_mean_N1e2,  res_norm_mean_N1e3,    res_norm_mean_N1e4, res_norm_mean_N1e5,     res_norm_mean_N1e6,    res_norm_mean_N1e7])
    
    
    points_o1var = [[1e-3, 1e-3], [1e-3, 5e-3], [5e-3, 5e-3]]   #label if  plotting as a function of 1/sqrt(N)
    points_o1var = [[1e-6, 1e-3], [1e-6, 1e-2], [1e-4, 1e-2]] 
    order= r'$\mathcal{O}(1/\sqrt{N})$'
    triangle_o1var = plt.Polygon(points_o1var, fill=None  ,edgecolor='grey') 
    
    Nlist = np.array([1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
    N_inv = 1/ Nlist
    plt.plot(N_inv,  res_norm_N, 'bo-')
    plt.xscale('log')
    plt.yscale('log')
   # plt.title('Tol = 1e-5')
    plt.gca().add_patch(triangle_o1var)
    plt.ylabel(r'$||  \hat{\rho^*} - \Phi^N_T(\hat{\rho^*)}||$'  , fontsize = 15)
    plt.xlabel(r'$1/N$', fontsize = 14)
    plt.annotate(order,  xy=(1.6e-6, 5e-3), xytext=(1.6e-6, 5e-3), fontsize=11, color='grey')
    plt.savefig("plots/Tolerance_on_NK-solution_converges_N-1.pdf")
    plt.show()
    
    
