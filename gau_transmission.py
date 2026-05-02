import numpy as np
from parameters import *
from density_solver import *


def gau_co_T(Ω0p, Ω0c, E0p, E0c, weights):

    """Computes transmission for co-propagating lasers with 
    inhomogeneous broadening in the EIT regime"""

    Ωp = Ω0p
    Ep = E0p
    Ωc = Ω0c
    Ec = E0c

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

        rp = np.exp(-αp* d / 2)
        rc = np.exp(-αc * d / 2)

        Ep *= rp
        Ec *= rc
        Ωp *= rp
        Ωc *= rc

    Tp = (Ep / E0p) ** 2
    Tc = (Ec / E0c) ** 2

    return Tp, Tc




def gau_counter_prop(Ω0p, E0p, Ω0c_guess, E0c_guess, weights):

    """Computes propagation for counter-propagating lasers with 
    inhomogeneous broadening in the EIT regime"""

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

    return Ωp, Ep, Ωc, Ec

def gau_counter_T(Ω0p, Ω0c, E0p, E0c, weights):

    """Computes transmission for counter-propagating lasers with 
    inhomogeneous broadening in the EIT regime using estimation of final 
    parameters of control laser"""

    Ec_guess = E0c
    Ωc_guess = Ω0c

    Ωp,Ep,Ωc,Ec = gau_counter_prop(Ω0p,E0p,Ωc_guess,Ec_guess, weights)

    while abs(Ec-E0c)/E0c>1e-7:

        # If the program takes too long to execute, invesitgate behavior of shooting method
        # print(abs(Ec-E0c)/E0c)
        
        Ec_guess += (E0c-Ec) 
        Ωc_guess = -dip2*Ec_guess/hbar

        Ωp,Ep,Ωc,Ec = gau_counter_prop(Ω0p,E0p,Ωc_guess,Ec_guess, weights)

    Tp = (Ep/E0p)**2 
    Tc = (Ec_guess/E0c)**2

    return Tp, Tc
