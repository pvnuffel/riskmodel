import numpy as np
import matplotlib.pylab as plt
import sys, glob

from tempfile import TemporaryFile
outfile = TemporaryFile()

#np.save(outfile, diff_norm )
#np.savetxt('results/data/diff_norm_x005.out', diff_norm)
#np.savetxt('results/data/Jv_sde_norm_x005.out', Jv_norm)

#np.savetxt('results/data/rho_x005.out', diff_norm)
#np.savetxt('results/data/Jv.out', Jv_norm)
Jv_norm = np.loadtxt('results/data/Jv.out')


#Jv_var = scipy.zeros((M))
#
#average= Jv_norm[-1]
#for m in range(1, len( Jv_norm)):
#    Jv_var= np.power((Jv_norm[m] - average),2) 
#
#plot_var = plt.plot(Jv_var)
#        
#np.savetxt('results/data/rho_norm_x001.out', rho_norm) 
rho42_005 = np.loadtxt('results/data/rho_42_x005.out')
rho42_001 = np.loadtxt('results/data/rho_42_x001.out')


datafile1= 'results/data/diff_norm_x001.out'
data1 = np.loadtxt(datafile1)

data1_systemic = data1[2500:]
system_error1= np.mean(data1_systemic)
print "system error for dx=0.01 = ", system_error1 

datafile2= 'results/data/diff_norm_x005.out'
data2 = np.loadtxt(datafile2)

data2_systemic = data2[2500:]
system_error2= np.mean(data2_systemic)
print "system error for dx=0.05 = ", system_error2 


ratio = np.power((0.01/0.05),-2)
print "Compare expected ratio = ", ratio ,  "with computed ratio = " , system_error1/system_error2   


rho_systemic= rho42_001[200:]
system_error1= np.mean(rho_systemic)
print "system error for dx=0.01 = ", system_error1 

datafile2= 'results/data/diff_norm_x005.out'
data2 = np.loadtxt(datafile2)

rho2_systemic = rho42_005[200:]
system_error2= np.mean(rho2_systemic)
print "system error for dx=0.05 = ", system_error2 


ratio = np.power((0.01/0.05),-2)
print "Compare expected ratio = ", ratio ,  "with computed ratio = " , system_error1/system_error2  


plot1 =plt.plot(rho_42_005, color='green', label = r'$\Delta x=0.05$')
plot2 =plt.plot(rho_42_001,  color ='blue', label = '$\Delta x=0.01$')
plt.legend(bbox_to_anchor=(1, 1), numpoints = 1 )

plt.xlabel('m', fontsize = 16)
plt.ylabel(r'$\Delta \rho$', fontsize = 20)

plt.annotate(r'$ \frac{\mathbb{E} \left[ \rho_{(\Delta x=0.01)} \right] }{\mathbb{E} \left[ \rho_{(\Delta x=0.05)} \right] }= 1.2028$', xy=(800, 0), xytext=(200, -0.1), fontsize=20)
            

plt.savefig("results/rho_particle42_ratio1p2.pdf")
plt.show()
plot1 =plt.plot(data1, color='green', label = r'$\Delta x=0.01$')
plot2 =plt.plot(data2,  label = '$\Delta x=0.05$')
#plot1 =plt.plot(np.append(np.roll(value,1),value[9]) )
#plot1 =plt.plot( norm_Jv_fp )
#M=10
#norm_Jv_fp 
#average_Jv_norm = scipy.zeros((M))
#for m in range(len(Jv_norm)):
    
    


ax=plt.subplot(111)
ax.set_xlim(1, 5000)


# if len(sys.argv) >2:
#       datafile2 = sys.argv[2]
#       datafile3 = sys.argv[3]
#      # datafile4 = sys.argv[4]
# else:
#       print "No data file name given. Please enter"
#       datafile = raw_input("-> ")
# if len(glob.glob(datafile))==0:
#       print "Data file %s not found. Exiting" % datafile
#       sys.exit() 
# data = np.loadtxt(datafile)
# data2 = np.loadtxt(datafile2)
# data3 = np.loadtxt(datafile3)

# line1 = pl.plot(data[:,0], data[:,1] , 'r',marker='o', markersize=3, label=r'$\bar{\nu}=0.5$')
# line2 = pl.plot(data2[:,0], data2[:,1], 'b',marker='o',  markersize=3, label=r'$\bar{\nu}=0.05$'  )
# line3 = pl.plot(data3[:,0], data3[:,1], 'g',marker='o',markersize=3,  label= r'$\bar{\nu}=0.5$'  )


#xmin =0
#xmax = 100
#ymin = 0
#ymax = 1
#plt.set_xlim(1,M)
plt.xscale('log')
plt.yscale('log')
#pl.ylabel(r'$T^*$', fontsize = 20)
plt.xlabel('m', fontsize = 16)
#plt.ylabel(r'$||\mathbb{E}(\mathbf{J}(\mathbf{U}_{SDE}) \cdot \mathbf{V})||_2 - ||\mathbf{J}(\mathbf{U}_{FP}) \cdot \mathbf{V})||_2 $', fontsize = 20)
plt.ylabel(r'$||\mathbb{E}\left[\mathbf{J}(\mathbf{U}_{SDE}) \cdot \mathbf{V} \right]- \mathbf{J}(\mathbf{U}_{FP}) \cdot \mathbf{V}  )||_2 $', fontsize = 20)
plt.legend([plot1, plot2], loc='best') #, ['x=0.05', 'x=0.01']) #, loc='best') #, numpoints=1)
plt.legend([plot1, plot2], loc='best') #, ['x=0.05', 'x=0.01']) #, loc='best') #, numpoints=1)
#pl.legend([line1, line2, line3],  'best')
#pl.legend( loc='best', numpoints = 1 )
#first_legend = pl.legend([line1], loc=1)
plt.legend(bbox_to_anchor=(1, 1), numpoints = 1 )
#pl.legend([line3], bbox_to_anchor=(0.3, 0.6))
#ax = pl.gca().add_artist(first_legend)
#plt.savefig("results/Jve_sde-Jv_pde_compare.pdf")



plt.show()
