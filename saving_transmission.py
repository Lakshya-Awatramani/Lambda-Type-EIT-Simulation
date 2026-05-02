from parameters import *
from gau_transmission import *
from tqdm import tqdm

# Power ratios are generated in log10 space

def T_varPc():

    """Saving probe and control transmissions, varying the control power"""

    filename = "var_Pc"

    Pratios = np.logspace(-5,0.1,ns) 


    Tpco = []; Tpcn = [] 
    Tcco = []; Tccn = [] 

    try:
        with open('npy/Tpco%s.npy' % filename, 'rb') as f: Tpco=np.load(f)
        with open('npy/Tpcn%s.npy' % filename, 'rb') as f: Tpcn=np.load(f)
        with open('npy/Tcco%s.npy' % filename, 'rb') as f: Tcco=np.load(f)
        with open('npy/Tccn%s.npy' % filename, 'rb') as f: Tccn=np.load(f)
    except (FileNotFoundError,NameError) as e:
        
        for rat in tqdm(Pratios):
            
            Pci = Pp/rat
            E0ci = (2*Pci/(c*ε0*Ac*n))**0.5
            Ω0ci = -dip2*E0ci/hbar

            temp = gau_co_T(Ω0p,Ω0ci,E0p,E0ci,weights)
            Tpco.append(temp[0])
            Tcco.append(temp[1])
            temp = gau_counter_T(Ω0p,Ω0ci,E0p,E0ci,weights)
            Tpcn.append(temp[0])
            Tccn.append(temp[1]) 

        with open('npy/Tpco%s.npy' % filename, 'wb') as f: np.save(f,Tpco)
        with open('npy/Tpcn%s.npy' % filename, 'wb') as f: np.save(f,Tpcn)
        with open('npy/Tcco%s.npy' % filename, 'wb') as f: np.save(f,Tcco)
        with open('npy/Tccn%s.npy' % filename, 'wb') as f: np.save(f,Tccn)


def T_varPp():

    """Saving probe and control transmissions, varying the probe power"""

    filename = "var_Pp"


    Pratios = np.logspace(-5,4,ns) 


    Tpco = []; Tpcn = [] 
    Tcco = []; Tccn = [] 
 
    try:
        with open('npy/Tpco%s.npy' % filename, 'rb') as f: Tpco=np.load(f)
        with open('npy/Tpcn%s.npy' % filename, 'rb') as f: Tpcn=np.load(f)
        with open('npy/Tcco%s.npy' % filename, 'rb') as f: Tcco=np.load(f)
        with open('npy/Tccn%s.npy' % filename, 'rb') as f: Tccn=np.load(f)
    except (FileNotFoundError,NameError) as e:
        

        for rat in tqdm(Pratios):
            
            Ppi = Pc*rat
            E0pi = (2*Ppi/(c*ε0*Ap*n))**0.5 
            Ω0pi = -dip1*E0pi/hbar 

            temp = gau_co_T(Ω0pi,Ω0c,E0pi,E0c,weights)
            Tpco.append(temp[0])
            Tcco.append(temp[1])
            temp = gau_counter_T(Ω0pi,Ω0c,E0pi,E0c,weights)
            Tpcn.append(temp[0])
            Tccn.append(temp[1]) 

        with open('npy/Tpco%s.npy' % filename, 'wb') as f: np.save(f,Tpco)
        with open('npy/Tpcn%s.npy' % filename, 'wb') as f: np.save(f,Tpcn)
        with open('npy/Tcco%s.npy' % filename, 'wb') as f: np.save(f,Tcco)
        with open('npy/Tccn%s.npy' % filename, 'wb') as f: np.save(f,Tccn)


T_varPc()
T_varPp()