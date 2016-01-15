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
    eps_list_exponents=[4, 5,6,7]
    eps_list= [1e-4, 1e-5, 1e-6, 1e-7]
    rho_Dt_fp= np.loadtxt('data/25-11-Rho_pde_dx_e-4_eps-5.out')   #,  np.loadtxt('25-11-Rho_pde_dx_e-4_eps-5.out') ,  
    #np.loadtxt('25-11-Rho_pde_dx_e-4_eps-6.out') ,  np.loadtxt('25-11-Rho_pde_dx_e-4_eps-7.out') ]  #these data should be the same
    Jv_pde= [ np.loadtxt('data/25-11-Jv_pde_dx_e-4_eps-4.out') ,  np.loadtxt('data/25-11-Jv_pde_dx_e-4_eps-5.out') ,  
    np.loadtxt('data/25-11-Jv_pde_dx_e-4_eps-6.out') ,  np.loadtxt('data/25-11-Jv_pde_dx_e-4_eps-7.out') ]
  #  Jv_pde = Jv_pde_list[1]
    
    #SDE       
   # sde_rho_sq = np.loadtxt('data/25-11-rho_sq_eps-5.out')
   # E_rho = np.loadtxt('data/25-11-E_rho_eps-5.out')
   # sq_E_rho  =np.loadtxt('data/25-11-sq_E_rho_eps-5.out')
    
    sde_Jv_sq = [ np.loadtxt('data/eps_fin-Jv_sq_eps-4.out'), np.loadtxt('data/eps_fin-Jv_sq_eps-5.out'), 
                         np.loadtxt('data/eps_fin-Jv_sq_eps-6.out'), np.loadtxt('data/eps_fin-Jv_sq_eps-7.out') ] 
    E_Jv = [np.loadtxt('data/eps_fin-E_Jv_eps-4.out'),  np.loadtxt('data/eps_fin-E_Jv_eps-5.out'), 
            np.loadtxt('data/eps_fin-E_Jv_eps-6.out'), np.loadtxt('data/eps_fin-E_Jv_eps-7.out') ]    # [eps_i x N_i x dx] = [4 x 7 x 3400]
    sq_E_Jv =  [np.loadtxt('data/eps_fin-sq_E_Jv_eps-4.out' ), np.loadtxt('data/eps_fin-sq_E_Jv_eps-5.out' ), 
                np.loadtxt('data/eps_fin-sq_E_Jv_eps-6.out' ),np.loadtxt('data/eps_fin-sq_E_Jv_eps-7.out' )]
                
 
    #rho_Dt_sde = E_rho[-1]
    Jv_sde = E_Jv[-1][-1]  #dim n_x

    
    resize_factor = int (len(rho_Dt_fp)/len(Jv_sde))
  #  resize_factor = 10
    print "Discretisation for solving sde is ",  resize_factor , " times coarser than the discretisation for solving the pde"
    
        
    def resize( original_vector, new_size):
        resize_factor = int (len(original_vector)/new_size) 
        print "Discretisation for solving sde is ",  resize_factor , " times coarser than the discretisation for solving the pde"  
        new_vector = scipy.zeros(new_size)
        for i in range (0,new_size):
            bin_av = 0
            for j in range (0,  resize_factor ):
                bin_av = bin_av+ original_vector[i*resize_factor+j]
                new_vector[i] = bin_av/resize_factor
        return new_vector

    
#    rho_coarse = scipy.zeros(len(rho_Dt_sde))
#    for i in range (0,len(rho_coarse)):
#        bin_av = 0
#        for j in range (0,  resize_factor ):
#            bin_av = bin_av+ rho_Dt_fp[i*resize_factor+j]
#        rho_coarse[i] = bin_av/resize_factor
            

    bins = len(Jv_sde)
    
    Jv_coarse = scipy.zeros((len(eps_list),bins))
    for eps_i in range (0,len(eps_list)):
        for i in range (0,len(Jv_coarse[eps_i])):
            bin_av = 0
            for j in range (0,  resize_factor ):
                bin_av = bin_av+ Jv_pde[eps_i][i*resize_factor+j]
                Jv_coarse[eps_i][i] = bin_av/resize_factor
                

                
                
#    pde_rho_sq = norm(rho_coarse)**2
                    
    pde_Jv_sq = scipy.zeros(len(eps_list))
    for eps_i in range (0,len(eps_list)):
          pde_Jv_sq[eps_i] =  norm(Jv_coarse[eps_i])**2
 
 
    Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000, 128000, 256000])
    Nlist_inv = 1.0/Nlist
    
    log_flag = 'log'
  #  log_flag = 'linear'
    save_flag = True
    
   # Nlist = scipy.array([1000,2000,4000,8000,16000])
       
       
#    plot_var = plt.plot(Nlist, (sde_rho_sq - sq_E_rho)/bins, label = '$Variance$')
#    plt.ylabel(r'Var$(\rho)$' , fontsize = 16)
#    #plt.ylabel()
#    plt.xlabel('$N$', fontsize = 16)
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    plt.show()

#    plot_mse = plt.plot(Nlist, (sde_rho_sq -pde_rho_sq)/bins, label = '$RMS^2$')
#    plt.xlabel('$N$', fontsize = 16)
#    plt.ylabel(r'MSE( $\rho$)' , fontsize = 16)
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    plt.show()
# 
#    plot_bias =  plt.plot(Nlist, (abs(pde_rho_sq - sq_E_rho))/bins, label = '$Bias$')
#    plt.xlabel('$N$', fontsize = 16)
#    plt.ylabel(r'Bias $(\rho)$', fontsize = 16)
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    plt.show()
    
 #   plot_rho = plt.plot( rho_coarse, label = '$\$')
 #   plot_rho =  plt.plot(rho_Dt_sde, label = '$\$')
  #  plt.ylabel(r' $\rho$', fontsize = 16)
 #   plt.xlabel('$x$', fontsize = 16)
  #  plt.show()

    lines = ["-","--","-.",":"]    
    
    points_o1var = [[1e-4, 2e-2], [1e-4, 1e-1], [5e-4, 1e-1]] 
 #   points_o1var = [[1e-5, 1e-2], [2e-5, 2e-2], [1e-5, 2e-2]] 
    triangle_o1var = plt.Polygon(points_o1var, fill=None  ,edgecolor='grey') 
    
    
    #VariancePlot
       
    for eps_i in range (0,len(eps_list)):
        plot_var = plt.plot(Nlist_inv, (sde_Jv_sq[eps_i] -sq_E_Jv[eps_i])/bins, lines[eps_i],
                            label =r'$\varepsilon=10^{-%d}$' %eps_list_exponents[eps_i] , linewidth=2)
    plt.ylabel('Var ($\mathbf{\hat{Jv}} $)', fontsize = 16)
    plt.xlabel('$1/N$', fontsize = 16) 
    plt.xscale(log_flag)
    plt.yscale(log_flag)
    plt.legend([plot_var], loc='best')
    plt.legend(bbox_to_anchor=(1, 0.5), numpoints = 1 )
    plt.gca().add_patch(triangle_o1var)
    if(save_flag): plt.savefig("plots/Var_N_eps.pdf")
    plt.show()
    
     
    #VBiasPlot
     
    bias =  scipy.zeros(((len(eps_list), len(Nlist), bins)))
    biasnorm=  scipy.zeros((len(eps_list),len(Nlist)))
    bias_sq=  scipy.zeros((len(eps_list), len(Nlist)))
    
    for eps_i in range (0,len(eps_list)):
        for i in range (0,len(Nlist)):
            bias[eps_i][i] = Jv_coarse[eps_i] - E_Jv[eps_i][i]
            biasnorm[eps_i][i] = norm(bias[eps_i][i])
            bias_sq[eps_i][i] = biasnorm[eps_i][i]**2
            
    for eps_i in range (0,len(eps_list)):
        plot_bias_sq = plt.plot(Nlist_inv, bias_sq[eps_i] /bins, lines[eps_i], label =r'$\varepsilon=10^{-%d}$' %eps_list_exponents[eps_i] )
    plt.xlabel('$1/N$', fontsize = 16)
    plt.ylabel(r'$\left(\mathbf{Bias}(\mathbf{\hat{Jv}} , \mathbf{\bar{Jv} }) \right)^2$'  , fontsize = 16)
    plt.xscale(log_flag)
    plt.yscale(log_flag)
    plt.legend([plot_bias_sq], loc='best')
    plt.legend(bbox_to_anchor=(1, 0.5), numpoints = 1 )
    if(save_flag): plt.savefig("plots/Bias_N_eps.pdf")
    plt.show()
    
    for eps_i in range (0,len(eps_list)):
        for i in range (0,len(Nlist)): 
            plot_bias= plt.plot( bias[eps_i][i], label=r'$N= %d$' %Nlist[i] )
        plt.xlabel('$n_x$', fontsize = 16)
        label =r'$\varepsilon=10^{-%d}$' %eps_list_exponents[eps_i] 
        plt.annotate(label,  xy=(60, 0.040), xytext=(1600, 0.6), fontsize=13)
        plt.ylabel( r'$\mathbf{Bias}(\mathbf{\hat{Jv}} , \mathbf{\bar{Jv} }) $', fontsize = 16)
        plt.legend([plot_bias], loc='best')
        plt.legend(bbox_to_anchor=(1, 1), numpoints = 1 )
        if(save_flag): plt.savefig('plots/Bias_eps__e-_%d.pdf' %eps_list_exponents[eps_i] )
        plt.show()

    #MSE Plot
    
    for eps_i in range (0,len(eps_list)):
        plot_Jv_mse = plt.plot(Nlist_inv, (sde_Jv_sq[eps_i] - 2*np.dot(E_Jv[eps_i],Jv_coarse[eps_i])  + pde_Jv_sq[eps_i])/bins, 
                               lines[eps_i], label =r'$\varepsilon=10^{-%d}$' %eps_list_exponents[eps_i] )
    plt.ylabel(r'MSE($\mathbf{\hat{Jv}} , \mathbf{\bar{Jv} }$)', fontsize = 16)
    plt.xlabel('$1/N$', fontsize = 16)
    plt.xscale(log_flag)
    plt.yscale(log_flag)
    plt.legend([plot_Jv_mse], loc='best')
    plt.legend(bbox_to_anchor=(1, 0.5), numpoints = 1 )
    if(save_flag): plt.savefig("plots/MSE_N_eps.pdf")
    plt.show()
    
    plot_Jv = plt.plot( Jv_coarse, label = '$\$')
    plot_Jv =  plt.plot(Jv_sde, label = '$\$')
    plt.ylabel(r' $Jv$', fontsize = 16)
   # plt.xlabel('$x$', fontsize = 16)
    plt.show()
    
    
                #plot_biasnorm = plt.plot(bias_sq[i], label=r'$N= %d$' %Nlist[i]  )
    
    
 
