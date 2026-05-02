import numpy as np
import matplotlib.pyplot as plt
from parameters import *

# Load files needed for plotting transmission

Tpco = np.load("npy/Tpcovar_Pc.npy")
Tpcn = np.load("npy/Tpcnvar_Pc.npy")
Tcco = np.load("npy/Tccovar_Pc.npy")
Tccn = np.load("npy/Tccnvar_Pc.npy")


Pratios = np.logspace(-5,0.1,ns) 

# PLotting the required graph for the loaded data

plt.ticklabel_format(scilimits=(0,0))
plt.plot(Pratios,Tpco,c='#90C9E3', linewidth=9,label=r'Probe: Co')
plt.plot(Pratios,Tpcn,c='#2481AC', linewidth=3,label=r'Probe: Counter')
plt.plot(Pratios,Tcco,c='#A4FF98',linewidth=9,label=r'Control: Co')
plt.plot(Pratios,Tccn,c='#00CF03',linewidth=3,label=r'Control: Counter')
plt.grid()
plt.legend(fontsize=15)
plt.title(r'Transmission',fontsize=20)
plt.xlabel(r'$P_p/P_c$',fontsize=20)
plt.xscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel(r'Transmission',fontsize=20)

plt.show()
