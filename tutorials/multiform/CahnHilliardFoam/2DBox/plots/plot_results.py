from plt_settings import *
import scipy.io as spio
#plt.style.use('seaborn')



# Import data from csv
gomez = np.loadtxt(open("GomezFine.csv", "rb"), delimiter=",", skiprows=0)
us = np.loadtxt(open("data_008.csv", "rb"), delimiter=",", skiprows=1)


plt.figure(figsize=(250 /25.4, 200 / 25.4))
plt.plot(us[:,4],us[:,0],'-r',label='OpenFOAM',linewidth=3)
plt.plot(gomez[:,0],gomez[:,1],'-b',label='Gomez et al. (2014)',linewidth=3)
#plt.plot([0.4999, 0.5001],[-1, 1], '--k', label='half way')
plt.ylabel('$\phi$',fontsize=30)
plt.xlabel('$x$',fontsize=30)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.ylim(-1, 1)
plt.xlim(0, 1)
#plt.show()
plt.legend(loc='best', fontsize=28)

#    plt.rc('grid', linestyle="-", color='black')
plt.grid(True)
#    plt.tight_layout()
plt.savefig('verification.pdf')
