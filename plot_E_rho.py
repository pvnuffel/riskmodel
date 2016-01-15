import numpy as np
import matplotlib.pylab as plt
import sys, glob

from tempfile import TemporaryFile
outfile = TemporaryFile()

#np.save(outfile, diff_norm )
#np.savetxt('results/data/rho_av_dx0p005', rho_av)
#np.savetxt('results/data/Jv_sde_norm_x005.out', Jv_norm)

#np.savetxt('results/data/rho_x005.out', diff_norm)
#np.savetxt('results/data/Jv.out', Jv_norm)
#Jv_norm = np.loadtxt('results/data/Jv.out')

#rho_av = np.loadtxt('results/data/rho_av_dx0p005')
#Jv_var = scipy.zeros((M))
#
#average= Jv_norm[-1]
#for m in range(1, len( Jv_norm)):
#    Jv_var= np.power((Jv_norm[m] - average),2) 
#
#plot_var = plt.plot(Jv_var)
#        
#np.savetxt('results/data/rho_norm_x001.out', rho_norm) 

plot1 =plt.plot( rho_norm/sqrt(len(rho_norm)))

ax=plt.subplot(111)
ax.set_xlim(1, 200)
plt.subplots_adjust(left=0.2) # move the left of the subplot to make room for the y-label

plt.annotate(r'$ \Delta x_{PDE} = 10^{-4}$ , $\Delta x_{SDE} = 5\cdot 10^{-3}$', xy=(60, 0.040), xytext=(42, 1.25), fontsize=13)
#plt.annotate(r'$ \Delta x_{PDE} = \Delta x_{SDE} = 0.005$', xy=(60, 0.15), xytext=(55, 0.20), fontsize=13)
plt.annotate(r'$ N=10000$', xy=(60, 0.15), xytext=(42, 1.1), fontsize=13)
#plt.xscale('log')
#plt.yscale('log')
#pl.ylabel(r'$T^*$', fontsize = 20)
plt.xlabel('$m$', fontsize = 18)
#plt.ylabel(r'$||\mathbb{E}(\mathbf{J}(\mathbf{U}_{SDE}) \cdot \mathbf{V})||_2 - ||\mathbf{J}(\mathbf{U}_{FP}) \cdot \mathbf{V})||_2 $', fontsize = 20)
#plt.ylabel(r'$||\mathbb{E}\left[\mathbf{J}(\mathbf{U}_{SDE}) \cdot \mathbf{V} \right]- \mathbf{J}(\mathbf{U}_{FP}) \cdot \mathbf{V}  )||_2 $', fontsize = 20)
plt.ylabel(r'$\frac{||\mathbb{E}[\mathbf{\rho}_{SDE}(m)  - \mathbf{\rho}_{FP}] ||_2 } {|| \mathbf{1}||_2}    $', fontsize = 18)

#plt.legend([plot1, plot2], loc='best') #, ['x=0.05', 'x=0.01']) #, loc='best') #, numpoints=1)
#pl.legend([line1, line2, line3],  'best')
#pl.legend( loc='best', numpoints = 1 )
#first_legend = pl.legend([line1], loc=1)
#plt.legend(bbox_to_anchor=(1, 1), numpoints = 1 )
#pl.legend([line3], bbox_to_anchor=(0.3, 0.6))
#ax = pl.gca().add_artist(first_legend)
plt.savefig("results/rho_norm_factor50.pdf")



plt.show()
