import matplotlib.pyplot as plt
import numpy as np
from density_solver import *
from tqdm import tqdm
from parameters import *


def gau_co_prop(Ω0p,Ω0c,E0p,E0c):

    """Saves laser parameters for co-propagation of lasers 
    in medium at EIT regime with inhomogeneous broadening"""

    Es = [] 

    Ωp = Ω0p
    Ep = E0p
    Ωc = Ω0c
    Ec = E0c 

    d = L/fidelity

    for _ in tqdm(range(fidelity)):

        Es.append([Ep,Ec]) 
        χp_eff = 0.0
        χc_eff = 0.0
        λp_eff = 0.0
        λc_eff = 0.0

        for Δ1_val, Δ2_val, w in zip(detunings, detunings, weights):

            λp = 2*pi*c/(k0p*c+Δ1_val) 
            λc = 2*pi*c/(k0c*c+Δ2_val)

            ρ=ρcl(Γ31=Γ31,Γ32=Γ32,γ2d=γ2d,γ3d=γ3d,Ωp=Ωp,Ωc=Ωc,Δ2=Δ2_val)
            ρsol = ρ.ρsol(Δ=Δ1_val) 
            ρ31=ρsol[6]; ρ32=ρsol[7]
            χp = N * dip1 * ρ31 / (ε0 * Ep)
            χc = N * dip2 * ρ32 / (ε0 * Ec)

            χp_eff += w * χp
            χc_eff += w * χc
            λp_eff += w * λp
            λc_eff += w * λc

        Ɛp1 = np.real(χp_eff) + 1
        Ɛp2 = np.imag(χp_eff)
        Ɛc1 = np.real(χc_eff) + 1
        Ɛc2 = np.imag(χc_eff)

        κp = Ɛp2/(np.sqrt(2*np.hypot(Ɛp1,Ɛp2)+Ɛp1))
        κc = Ɛc2/(np.sqrt(2*np.hypot(Ɛc1,Ɛc2)+Ɛc1))

        αp = (4*pi*κp)/λp_eff 
        αc = (4*pi*κc)/λc_eff

        rp = np.exp(-αp* d / 2)
        rc = np.exp(-αc * d / 2)

        Ep *= rp
        Ec *= rc
        Ωp *= rp
        Ωc *= rc

    return np.array(Es) 




def gau_counter_prop(Ω0p,Ω0c_guess,E0p,E0c_guess):

    """Saves laser parameters for counter-propagation of lasers 
    in medium at EIT regime with inhomogeneous broadening"""


    Ωp = Ω0p
    Ep = E0p
    Ωc = Ω0c_guess
    Ec = E0c_guess

    d = L / fidelity 


    for _ in range(fidelity):

        χp_eff = 0.0
        χc_eff = 0.0
        λp_eff = 0.0
        λc_eff = 0.0

        for Δ1_val, Δ2_val, w in zip(detunings, detunings, weights):
            
            λp = 2*pi*c/(k0p*c+Δ1_val) 
            λc = 2*pi*c/(k0c*c+Δ2_val)

            ρ=ρcl(Γ31=Γ31,Γ32=Γ32,γ2d=γ2d,γ3d=γ3d,Ωp=Ωp,Ωc=Ωc,Δ2=Δ2_val)
            ρsol = ρ.ρsol(Δ=Δ1_val) 
            ρ31=ρsol[6]; ρ32=ρsol[7]

            χp = N * dip1 * ρ31 / (ε0 * Ep)
            χc = N * dip2 * ρ32 / (ε0 * Ec)

            χp_eff += w * χp
            χc_eff += w * χc
            λp_eff += w * λp
            λc_eff += w * λc

        Ɛp1 = np.real(χp_eff) + 1
        Ɛp2 = np.imag(χp_eff)
        Ɛc1 = np.real(χc_eff) + 1
        Ɛc2 = np.imag(χc_eff)

        κp = Ɛp2/(np.sqrt(2*np.hypot(Ɛp1,Ɛp2)+Ɛp1))
        κc = Ɛc2/(np.sqrt(2*np.hypot(Ɛc1,Ɛc2)+Ɛc1))

        αp = (4*pi*κp)/λp_eff 
        αc = (4*pi*κc)/λc_eff

        rp = np.exp(-αp * d / 2)
        rc = np.exp(αc * d / 2)

        Ep *= rp
        Ec *= rc
        Ωp *= rp
        Ωc *= rc

    return Ωc,Ec
    

def gau_counter_prop_save(Ω0p,Ω0c,E0p,E0c):

    """Computes propagation for counter-propagating lasers with inhomogeneous broadening 
    in the EIT regime using estimation of final 
    parameters of control laser"""

    Ec_guess = E0c
    Ωc_guess = Ω0c


    Ωc,Ec = gau_counter_prop(Ω0p,Ωc_guess,E0p,Ec_guess)

    while abs(Ec-E0c)/E0c>1e-7:
        
        # If the program takes too long to execute, invesitgate behavior of shooting method
        print(abs(Ec-E0c)/E0c)

        Ec_guess += (E0c-Ec) 
        Ωc_guess = -dip2*Ec_guess/hbar 

        Ωc,Ec = gau_counter_prop(Ω0p,Ωc_guess,E0p,Ec_guess)
    

    Es = []
    Ωp = Ω0p
    Ep = E0p
    Ωc = Ωc_guess
    Ec = Ec_guess

    d = L / fidelity 


    for _ in tqdm(range(fidelity)):

        χp_eff = 0.0
        χc_eff = 0.0
        λp_eff = 0.0
        λc_eff = 0.0

        for Δ1_val, Δ2_val, w in zip(detunings, detunings, weights):
            
            λp = 2*pi*c/(k0p*c+Δ1_val) 
            λc = 2*pi*c/(k0c*c+Δ2_val)

            ρ=ρcl(Γ31=Γ31,Γ32=Γ32,γ2d=γ2d,γ3d=γ3d,Ωp=Ωp,Ωc=Ωc,Δ2=Δ2_val)
            ρsol = ρ.ρsol(Δ=Δ1_val) 
            ρ31=ρsol[6]; ρ32=ρsol[7]

            χp = N * dip1 * ρ31 / (ε0 * Ep)
            χc = N * dip2 * ρ32 / (ε0 * Ec)

            χp_eff += w * χp
            χc_eff += w * χc
            λp_eff += w * λp
            λc_eff += w * λc

        Ɛp1 = np.real(χp_eff) + 1
        Ɛp2 = np.imag(χp_eff)
        Ɛc1 = np.real(χc_eff) + 1
        Ɛc2 = np.imag(χc_eff)

        κp = Ɛp2/(np.sqrt(2*np.hypot(Ɛp1,Ɛp2)+Ɛp1))
        κc = Ɛc2/(np.sqrt(2*np.hypot(Ɛc1,Ɛc2)+Ɛc1))

        αp = (4*pi*κp)/λp_eff 
        αc = (4*pi*κc)/λc_eff

        rp = np.exp(-αp * d / 2)
        rc = np.exp(αc * d / 2)

        Ep *= rp
        Ec *= rc
        Ωp *= rp
        Ωc *= rc

        Es.append([Ep, Ec])


    return np.array(Es) 



Esco = gau_co_prop(Ω0p,Ω0c,E0p,E0c)
Escn = gau_counter_prop_save(Ω0p,Ω0c,E0p,E0c)


x = np.linspace(0,L,fidelity)
fig, (ax1, ax2) = plt.subplots(1, 2)

fig.set_figheight(5)
fig.set_figwidth(13)
fig.set_dpi(30)

ax1.ticklabel_format(scilimits=(0,0))
ax2.ticklabel_format(scilimits=(0,0))

ax1.plot(x,Esco[:,0],c='orangered',linewidth=7,label=r'Co')
ax1.plot(x,Escn[:,0],c='steelblue',linewidth=3,label=r'Counter')
ax1.set_xlabel(r'Distance along the cell (m)')
ax1.set_ylabel(r'$E_p$ (V/m)')
ax1.set_title(r'Probe')
ax1.legend()

ax2.plot(x,Esco[:,1],c='orangered',linewidth=3,label=r'Co')
ax2.plot(x,Escn[:,1],c='steelblue',linewidth=3,label=r'Counter')
ax2.set_xlabel(r'Distance along the cell (m)')
ax2.set_ylabel(r'$E_c$ (V/m)')
ax2.set_title(r'Control')
ax2.legend()

ax1.grid()
ax2.grid()

fig.suptitle(r'Propagation: $L=%.1f$ cm, $P_p=%.1f$μW, $P_c=%.1f$ mW' % (L*100,Pp*1e6,Pc*1e3))

plt.show()
