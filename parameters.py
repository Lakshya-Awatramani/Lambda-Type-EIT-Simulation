from scipy import constants
import numpy as np
from scipy.stats import norm


# Constants in SI units
pi = constants.pi
ε0 = constants.epsilon_0
hbar = constants.hbar 
c = 299792458 

fidelity = int(1e3)

# Control & probe laser diameter
Dc = 2e-3 
Dp = 2e-3

# Laser spot size of control & probe (m^2)
Ac = pi*(Dc/2)**2 
Ap = pi*(Dp/2)**2 

# Average index of refraction
n = 1

# Control & probe laser power (W)
Pc = 1e-3
Pp = 20.0e-6

# Amplitude of control & probe field (V/m)
E0c = (2*Pc/(c*ε0*Ac*n))**0.5 
E0p = (2*Pp/(c*ε0*Ap*n))**0.5

# D1 transition dipole moment (Cm)
dip = -2.537e-29  
dip1 = (1/4)**0.5*dip 
dip2 = (1/12)**0.5*dip 

# Vacuum wavelength (m) and wavenumber (1/m) of control & probe field
λ0c = 795.0e-9
λ0p = 795.0e-9
k0c = 2*pi/λ0c
k0p = 2*pi/λ0p

# Control & probe Rabi frequency (rad Hz)
Ω0c = -dip2*E0c/hbar
Ω0p = -dip1*E0p/hbar


# Spontaneous emission rate out of state |3>, |3>->|2> & |3>->|1> (rad Hz) 
Γ3 = 36.1e6 * 2 *np.pi 
Γ32 = Γ3/2 
Γ31 = Γ3/2 

# Dephasing rate of state |3> & |2> (Hz)
γ3d = 0.1 * Γ3 
γ2d = 1e2 * 2 *np.pi


# Detuning of probe & control field (Hz)
Δ1 = 0 
Δ2 = 0 

# Atomic number density of the Rb-87 sample (1/m^3)
N = 9.715e15 
# Length cell (m)
L = 5e-2
# Power ratio distribution factor
ns = 10


# Inhomogeneous broadening setup (gaussian profile)

detunings = np.linspace(-100e6, 100e6, 41) * 2*np.pi
mu = 0
sigma = 5e6 * 2*np.pi
weights = norm.pdf(detunings, mu, sigma)
weights /= np.sum(weights)
