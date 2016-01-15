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
from matplotlib.path import Path #to draw trriangles indicating order of convergence
import matplotlib.patches as patches #to draw trriangles indicating order of convergence

from scipy.linalg import norm
from scipy import sqrt

import matplotlib.pylab as plt


if __name__=="__main__":
   
   
   
   

    Mlist= [ 100, 200, 300, 400, 500]
    
   # Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000]) #,128000])
  #  Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000, 128000, 256000])
    Nlist=scipy.array([1000])
    Nlist_inv = 1.0/Nlist
    
        #SDE       

    E_rho = scipy.zeros(((len(Mlist), 34)))    # [m x N_i x dx] = [5 x 8 x 3400]
    for M_factor in range(0,len(Mlist)):
        m= (M_factor+1)*100
        E_rho[M_factor] = np.loadtxt('test-E_rho_m%d.out' %m)
#    sq_E_rho = scipy.zeros((len(Mlist), len(Nlist))) 
#    for M_factor in range(0,5):
#        m= (M_factor+1)*20
#        sq_E_rho[M_factor] = np.loadtxt('data/8-12-sq_E_rho_m%d.out' %m)

#    sde_Jv_sq =scipy.zeros((len(Mlist), len(Nlist)))
#    for M_factor in range(0,len(Mlist)):
#        m= (M_factor+1)*20
#        sde_Jv_sq[M_factor] = np.loadtxt('data/fin-Jv_sq_m%d.out' %m)
    E_Jv = scipy.zeros(((len(Mlist), 34)))    # [m x N_i x dx] = [5 x 8 x 3400]
    for M_factor in range(0,len(Mlist)):
        m= (M_factor+1)*100
        E_Jv[M_factor] = np.loadtxt('test-E_Jv_m%d.out' %m)     
#    sq_E_Jv = scipy.zeros((len(Mlist), len(Nlist)))  
#    for M_factor in range(0,len(Mlist)):
#        m= (M_factor+1)*20
#        sq_E_Jv[M_factor] = np.loadtxt('data/fin-sq_E_Jv_m%d.out' %m) 
        
#   rho_Dt_sde = E_rho[-1][-1]
    Jv_sde = E_Jv[-1]  #dim n_
    bins = len(Jv_sde)
    rho_sde = E_rho[-1]
    
    rho_sde = np.loadtxt('../../Particles/checkSystem/E_rho_M100_nw.out')
    Jv_sde= np.loadtxt('../../Particles/checkSystem/E_Jv_M100_nw.out')
    
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


    
    
    rho_Dt_fp= np.loadtxt('data/25-11-Rho_pde_dx_e-4_eps-5.out')   #,  np.loadtxt('25-11-Rho_pde_dx_e-4_eps-5.out') ,  
    #np.loadtxt('25-11-Rho_pde_dx_e-4_eps-6.out') ,  np.loadtxt('25-11-Rho_pde_dx_e-4_eps-7.out') ]  #these data should be the same
    Jv_pde= np.loadtxt('data/25-11-Jv_pde_dx_e-4_eps-5.out')     
  
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0+0.08, box.y0, box.width*0.95, box.height]) 
    #r = scipy.zeros(10)  
    dxlist= [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]
    bias_sq_dxlist = scipy.zeros(len(dxlist))    
    bias_sq_rho = scipy.zeros(len(dxlist))  
    for dx_i in range(0, len(dxlist)):
      #  bias_sq_rlist[r_inv-1] =  norm(resize( np.loadtxt('Jv_dx=e-2_r_inv=%d.out' %r_inv)  , bins) - Jv_sde )**2  Jv_dx=e-2_r=%d.out
        bias_sq_dxlist[dx_i] =  norm(resize( np.loadtxt('Jv_dx=%.3f.out' %dxlist[dx_i])  , bins) - Jv_sde )**2
        bias_sq_rho[dx_i] =  norm(resize( np.loadtxt('rho_dx=%.3f.out' %dxlist[dx_i])  , bins) -  rho_sde  )**2
      #  r[r_inv-1]= 1.0/r_inv
    plt.plot( dxlist, bias_sq_dxlist /bins, linewidth=2  ) 
    plt.xlabel('$\Delta x$', fontsize = 14)
    plt.ylabel(r'$\left(\hat{\mathbf{Bias}}(\mathbf{\hat{Jv}} , \mathbf{Jv}_{FP} \right))^2$'  , fontsize = 15)
   # plt.xscale(log_flag)
   # plt.yscale(log_flag)
   # plt.ticklabel_format(style='sci',  scilimits=(0,10))
  #  plt.ticklabel_format(style='sci') 
    plt.ticklabel_format(useOffset=False, axis='y') 
 #   plt.savefig('plots/ftcs_checkup.pdf')
    #plt.savefig('Jv_bias_dx_nw.pdf')
    plt.show()


    plt.plot(dxlist,  bias_sq_rho/bins, linewidth=2  )
    plt.xlabel('$\Delta x$', fontsize = 14)
    plt.ylabel(r'$\left(\hat{\mathbf{Bias}}(\mathbf{\hat{\rho}} , \mathbf{\rho}_{FP} \right))^2$'  , fontsize = 15)

   # plt.savefig('rho_bias_dx_nw.pdf')
    plt.show()
  #  Jv_pde = Jv_pde_list[1]


 
    xL = -1.7
    xR = 1.7
    dx=1e-3
    grid = scipy.arange(xL+dx/2.,xR,dx)
    
    v=scipy.zeros_like(grid)
    for j in range(len(grid)): 
        v[j]= np.sin(j*2*np.pi/len(grid))
        
    
 
    
    rho_coarse = resize(rho_Dt_fp,bins)
#    Jv_coarse_dxem2_r2 = resize( Jv_pde_dxem2_r0p001, bins )
    Jv_coarse = resize(Jv_pde,bins)
                
  #  pde_rho_sq = norm(rho_coarse)**2
    pde_Jv_sq = norm(Jv_coarse)**2

    
    log_flag = 'log'
  #  log_flag = 'linear'
    save_flag = False
    
    lines = ["-","--","-.",":","-","--","-.",":","-","--","-.",":"]    

   # points_o1a = [[1e-4, 1e-3], [1e-4, 5e-3], [5e-4, 5e-3]] 
    points_o1b = [[1e-5, 1e-3], [1e-5, 5e-3], [5e-5, 5e-3]] 
   # triangle_o1a = plt.Polygon(points_o1a, fill=None  ,edgecolor='r') 
    triangle_o1b = plt.Polygon(points_o1b, fill=None  ,edgecolor='grey')      

#Plot RHO
 #   plot_rho =  plt.plot(rho_Dt_sde, '-', color='g', label = '$\$')
#    plot_rho = plt.plot( rho_coarse, '--', color='r', label = '$\$')

   # plt.ylabel(r' $\rho$', fontsize = 16)
 #   plt.xlabel('$x$', fontsize = 16)
  #  plt.show()

        
    points_o1 = [[1e-5, 1e-3], [1e-5, 4e-3], [4e-5, 4e-3]]  
    triangle_o1 = plt.Polygon(points_o1, fill=None  ,edgecolor='r')
  
#
#    for m in range (0,len(Mlist)):
#        plot_var = plt.plot(Nlist_inv, (sde_rho_sq[m] -sq_E_rho[m])/bins, lines[m], label =r'$m=%d$' %Mlist[m] , linewidth=3)
#    plt.ylabel(r'Var ($\mathbf{\hat{\rho}} $)', fontsize = 14)
#    plt.xlabel('$1/N$', fontsize = 14) 
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    plt.legend([plot_var], loc='best')
#    plt.legend(bbox_to_anchor=(1, 0.5), numpoints = 1 )
#    if(save_flag): plt.savefig("plots/25-11_rho_Var_N_m.pdf")
#    plt.show()
#    
    #VBiasPlot
     
    bias =  scipy.zeros((len(Mlist),  bins))
    biasnorm=  scipy.zeros((len(Mlist)))
    bias_sq_rho=  scipy.zeros(len(Mlist))
    
    for m in range (0,len(Mlist)):
            bias[m] = rho_coarse - E_rho[m]
            biasnorm[m] = norm(bias[m])
            bias_sq_rho[m] = biasnorm[m]**2
            

    plot_bias_sq = plt.plot(Mlist, bias_sq_rho/bins, lines[1], label ='\Delta x = 0.1', linewidth=2  ) #
    plt.xlabel('$M$', fontsize = 16)
    plt.ylabel(r'$\left(\mathbf{Bias}(\mathbf{\hat{\rho}} , \mathbf{\bar{\rho} }) \right)^2$'  , fontsize = 16)
    
 #   plt.gca().add_patch(triangle_o1)
   # plt.gca().add_patch(triangle_o1)
    plt.annotate(r'\mathcall{O}(1)$', xy=(1e-5, 1e-3), xytext=(55, 0.125), fontsize=13)
   # plt.xscale(log_flag)
    plt.yscale(log_flag)
    plt.legend([plot_bias_sq], loc='best')
    plt.legend(bbox_to_anchor=(0.3, 1), numpoints = 1 )
    if(save_flag): plt.savefig("plots/25-11_rho_Bias_N.pdf")
    plt.show()
    
#PLOT Jv  
    
    #VariancePlot
    
#    points_o1var = [[1e-4, 2e-2], [1e-4, 1e-1], [5e-4, 1e-1]] 
# #   points_o1var = [[1e-5, 1e-2], [2e-5, 2e-2], [1e-5, 2e-2]] 
#    triangle_o1var = plt.Polygon(points_o1var, fill=None  ,edgecolor='grey') 
#    
#    
#    ax = plt.subplot(111)
#    box = ax.get_position()
#    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])   
#    for m in range (0,len(Mlist)):
#        plot_var = plt.plot(Nlist_inv, (sde_Jv_sq[m] -sq_E_Jv[m])/bins, lines[m], label =r'$M=%d$' %Mlist[m] , linewidth=2  )
#    plt.ylabel('Var ($\mathbf{\hat{Jv}} $)', fontsize = 16)
#    plt.xlabel('$1/N$', fontsize = 16) 
#    plt.gca().add_patch(triangle_o1var)
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    plt.legend([plot_var], loc='best')
#    plt.legend(bbox_to_anchor=(1, 0.5), numpoints = 1 )
#    if(save_flag): plt.savefig("plots/25-11Var_N_m.pdf")
#    plt.show()

    #VBiasPlot

    bias =  scipy.zeros(((len(Mlist),  bins)))
    biasnorm=  scipy.zeros((len(Mlist)))
    bias_sq = scipy.zeros(len(Mlist))
    
    for m in range (0,len(Mlist)):
        bias[m] = Jv_coarse - E_Jv[m]
        biasnorm[m] = norm(bias[m])
        bias_sq[m] = biasnorm[m]*2
            
    bias =  scipy.zeros(((len(Mlist),  bins)))
    biasnorm=  scipy.zeros((len(Mlist)))
    bias_sq_dxem1=  scipy.zeros(len(Mlist))
    
    


    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plot_bias_sq = plt.plot(Mlist, bias_sq /bins, lines[1], label ='$\Delta x = 10^{-4}$', linewidth=2  ) 

    plt.xlabel('$M$', fontsize = 16)
    plt.ylabel(r'$\left(\hat{\mathbf{Bias}}(\mathbf{\hat{Jv}} , \mathbf{Jv}_{FP} \right)^2$'  , fontsize = 16)
   # plt.xscale(log_flag)
    plt.yscale(log_flag)
    plt.legend([plot_bias_sq], loc='best')
   # plt.legend(bbox_to_anchor=(0.32, 1), numpoints = 1 )
    plt.legend(bbox_to_anchor=(1.43, 0.9), numpoints = 1 )
          
    if(save_flag): plt.savefig("plots/Bias_Jv_M_N.pdf")
  #  plt.savefig("plots/Bias_Jv_M_N.pdf")
    plt.show()
    
    
#    for m in range (0,len(Mlist)):
#        for i in range (0,len(Nlist)): 
#            plot_bias= plt.plot( bias[m][i], label=r'$N= %d$' %Nlist[i] )
#        plt.xlabel('$n_x$', fontsize = 16)
#        label =r'$m=%d$' %Mlist[m] 
#        plt.annotate(label,  xy=(20, 0.040), xytext=(160, 0.15), fontsize=13)
#        plt.ylabel( r'$\hat{\mathbf{Bias}}(\mathbf{\hat{Jv}} ,  \mathbf{Jv}_{FP} }) $', fontsize = 16)
#        plt.legend([plot_bias], loc='best')
#        plt.legend(bbox_to_anchor=(1, 1), numpoints = 1 )
#        if(save_flag): plt.savefig('plots/25-11Bias_e-_%d.pdf' %Mlist_exponents[m] )
#        plt.show()

    
 
