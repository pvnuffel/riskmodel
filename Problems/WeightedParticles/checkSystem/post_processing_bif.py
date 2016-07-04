# -*- coding: utf-8 -*-
import scipy
import numpy as np
from pylab import *
from scipy.linalg import norm



def get_fixpoint( mu, sigma, alpha):
    xi=0
    C=  sigma**2/(2*alpha)
    if (3*C>1):
        print "Make sure that 3D/(2alpha) <1 to get non-zero equilibrium solutions"
    else:
        xi = sqrt(1-3*C)*(1+ 6*mu/(sigma**2)*C**2*(1-2*C)/(1-3*C)) 
    return xi

def get_critical_alpha( mu, sigma):
    xi=0
    C=  sigma**2/(2*alpha)
    if (3*C>1):
        print "Make sure that 3D/(2alpha) <1 to get non-zero equilibrium solutions"
    else:
        xi = sqrt(1-3*C)*(1+ 6*mu/(sigma**2)*C**2*(1-2*C)/(1-3*C)) 
    return xi



if __name__=="__main__":
    
    dx=1e-2
    xL = -1.7
    xR = 1.7
    grid = scipy.arange(xL+dx/2.,xR,dx)
    
    sigma=1.
    mu = 0.1
    
    
        
    alphas = scipy.arange(0, 6, 0.11)
    fixpoints = scipy.zeros_like(alphas)
    for i in range (len(fixpoints)):
        fixpoints[i] = get_fixpoint(mu, sigma, alphas[i])
    
    plot_bif = plt.plot(alphas, fixpoints, 'b')
    plot_bif = plt.plot(alphas, -fixpoints, 'b')
    plt.ylabel(r'$\xi$', fontsize = 18)
    plt.xlabel(r'$\alpha$', fontsize = 16)
    plt.title(r'$\mu=0.1$ , $\sigma=1$')  
    plt.savefig('plots/bif_alpha.pdf')
    plt.show()
    
    
    #STOCHASTIC

    save_flag= True
    branch_N1e4 = np.loadtxt('Biff_Ne4_80steps.txt')
    branch_N1e5 = np.loadtxt('Biff_Ne5_80steps.txt')
    branch_N1e5_RL = np.loadtxt('Biff_Ne5_80steps_RL.txt')

    #branch_N1e5 = np.loadtxt('Biff_Ne5_tol_alpha_step_1_tol5e-3.txt')
    #branch_N1e5_l = np.loadtxt('Biff_Ne5_alpha_tol4e-3.txt')
  #  branch_N1e5_n = np.loadtxt('Biff_Ne5_alpha_tol3e-3.txt')


  #  branch_N1e6 = np.loadtxt('Biff_Ne6_tol_alpha_step_1_tole-3.txt')
    branch_N1e6_RL = np.loadtxt('Biff_Ne6_40steps_tol1e-3.txt')
    branch_N1e6_LR = np.loadtxt('Biff_Ne6_20steps__LR_tol1e-3.txt')
    branch_N1e6_LR = -np.loadtxt('Biff_Ne6_80steps.txt')
    branch_N1e6_RL = np.loadtxt('Biff_Ne5_80steps_RL.txt')
    

    sigma=1.
    mu = 0.1
    alpha= 2.

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
    
    alphas = scipy.arange(0, 8.3, 0.11)
    fixpoints = scipy.zeros_like(alphas)
    for i in range (len(fixpoints)):
        fixpoints[i] = get_fixpoint(mu, sigma, alphas[i])
    
    plot_bif = plt.plot(alphas, fixpoints, 'b', linewidth=1.7, label = r'$\xi$' )
    plot_bif = plt.plot(alphas, -fixpoints, 'b', linewidth=1.7)  #label = r'$\xi$' )
    plt.ylabel(r'$\bar{x}^*$', fontsize = 18)
    plt.xlabel(r'$\alpha$', fontsize = 16)
  #  plt.savefig('plots/bif_alpha.pdf')

    alpha_list= scipy.arange(2, 9, 1)  #voor LR
    alpha_list_l= scipy.arange(0, 3.5, 0.5)  #voor LR
    alpha_list_n = scipy.arange(8, -0.5, -0.5)  #voor RL
    alpha_list_LR = scipy.arange(0, 8.1, 0.4)  #voor LR
    alpha_list_RL = scipy.arange(8.1, 0, -0.2)  #voor RL
    alpha_list_80_LR = scipy.arange(0, 8.0001, 0.1)  #voor RL
    alpha_list_80_RL = scipy.arange(8.0001, 0, -0.1)  #voor RL
    
   # alpha_list =  scipy.arange(8, -0.5, -0.5)  #voor RL
    
  #  plt.plot(alpha_list, branch_N1e5, label=r'$N=10^5$')
  #  plt.plot(alpha_list_l, branch_N1e6, label=r'$N=10^6$')
  #  plt.plot(alpha_list_l, branch_N1e5_l)
   # plt.plot(alpha_list_n, branch_N1e5_n, 'r', label=r'$N=10^5$')
    
    plt.legend(bbox_to_anchor=(1, 0.75), numpoints = 1,  prop={'size':11})
   # plt.xlim([0,10])
    plt.savefig('plots/bifurcation_0.pdf')
    
    plt.plot(alpha_list_80_LR, -branch_N1e5 , 'm', label=r'$ N=10^5 $ ') 
    x = alpha_list_80_LR
    y= -branch_N1e5
    plt.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1], width=0.0043, scale_units='xy', angles='xy', scale=1, color='m')
    plt.plot(alpha_list_80_RL, -branch_N1e5_RL , 'y', label=r'$N=10^5$ ') 
    x = alpha_list_80_RL
    y= -branch_N1e5_RL
    plt.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1], width=0.0043, scale_units='xy', angles='xy', scale=1, color='y')
    
    plt.legend(bbox_to_anchor=(1, 0.75), numpoints = 1,  prop={'size':11})
    plt.savefig('plots/bifurcation_1.pdf')
    
    plt.plot(alpha_list_80_LR,  branch_N1e6_LR , 'g', label=r'$N=10^6$ ') 
    plt.quiver(alpha_list_80_LR[:-1],  branch_N1e6_LR[:-1] , alpha_list_80_LR[1:] -alpha_list_80_LR[:-1] ,  branch_N1e6_LR[1:]- branch_N1e6_LR[:-1] , 
               scale_units='xy', angles='xy', width=0.0043, scale=1, color='g')
    plt.plot(alpha_list_80_RL,  branch_N1e6_RL , 'r', label=r'$N=10^6$ ') 
    plt.quiver(alpha_list_80_RL[:-1],  branch_N1e6_RL[:-1] , alpha_list_80_RL[1:] -alpha_list_80_RL[:-1] ,  branch_N1e6_RL[1:]- branch_N1e6_RL[:-1] , 
               scale_units='xy', width=0.0043, angles='xy', scale=1, color='r')
               
 
             
   # plt.plot(alpha_list_80, branch_N1e4 , 'k', label=r'$   N=10^4$') 

   
    plt.legend(bbox_to_anchor=(1, 0.75), numpoints = 1,  prop={'size':11})
    plt.savefig('plots/bifurcation_fin.pdf')
    plt.show()
#    
#    x = alpha_list_80
#    
#    y =-branch_N1e5_RL 
#
#    plt.figure()
#    plt.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1], width=0.0043, scale_units='xy', angles='xy', scale=1)
#
#    plt.show()
    
    
    
##sigma_list= scipy.arange(2, 0.9, -0.1)   #voor RL
#sigma_list= scipy.arange(1.0, 2.1, 0.1)  #voor LR
#        
#for i in range (N_points):
#    if (i % 2 == 0): 
#        par_label = r'$\sigma= %g $' % sigma_list[i]
#        plt.plot(grid, rho_mean[i], label=par_label)   
#plt.xlim([-1.7,1.7])
#plt.ylim([0,0.6])
#plt.xlabel('$x$', fontsize = 12)
#plt.ylabel(r'$\rho^*$', fontsize = 12)
#plt.legend(bbox_to_anchor=(1,0.6),  numpoints = 1 )
#plt.legend(prop={'size':9})
#if save_flag:
#    plt.savefig('plots/bif/fixed_states_sde(sigma)_Ne6_mean_M10_LR.pdf')  
#plt.show()
#
##PLOT BIFURCATION 
#
#ax = plt.subplot(111)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.95, box.height])       
#plt.plot( sigma_list_fine,y_anal_fine,  'r', markersize=4,label= 'Analytic solution')
#plt.errorbar(sigma_list,  y_mean, yerr=y_error,   fmt='bo', ecolor='b', markersize=3, label='Newton-Krylov solution')
#plt.ylabel(r'$\left\||\rho^*\right\||_2$', fontsize = 12)
#plt.xlabel(r'$\sigma$', fontsize = 11)
#plt.xlim([0.9,2.1])
#plt.legend(prop={'size':8})
#
##plt.plot(x_sde, y_sde, 'go', markersize=4)
#
#if save_flag:
#    plt.savefig('plots/bif/bifurcation_sde_Ne6_anal(sigma)_LR.pdf')
#
#plt.show()

#
#plt.xlim([0.9,2.1])
#plt.ylabel(r'Bias($\left\||\rho^*\right\||_2)$', fontsize = 12)
#plt.xlabel(r'$\sigma$', fontsize = 12)
##plt.yscale('log')
#plt.errorbar(sigma_list, y_mean-y_anal, yerr=y_error,   fmt='bo', ecolor='b', markersize=3) # Mean of M realizations
#
#if save_flag:
#    plt.savefig('plots/bif/bifurcation_bias(sigma)_sde_Ne6_mean_M10_LR.pdf')
#
#plt.show()