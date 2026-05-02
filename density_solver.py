import numpy as np
from numpy import conjugate as conj
import matplotlib.pyplot as plt

# Density matrix element
ρijs = ['ρ11','ρ12','ρ13','ρ21','ρ22','ρ23','ρ31','ρ32','ρ33'] 


class ρcl:

    def ρsol(self,Δ):

        """Solves the overdetermined system to produce density matrix elements"""
        
        L3E2 = -1.j*(self.Δ2-Δ)+1/2*self.γ21
        L4E3 =  1.j*Δ+1/2*self.γ31
        L5E4 =  1.j*(self.Δ2-Δ)+1/2*self.γ21
        L7E6 =  1.j*self.Δ2+1/2*self.γ32
        L7E8 = -1.j*Δ+1/2*self.γ31
        L7E9 =  1.j*self.Δ2+1/2*self.γ32
        a = 1/2*1.j*self.Ωp; b = 1/2*1.j*self.Ωc

        # Overdetermined system
        ρdot = np.array([   [1,0,0,0,1,0,0,0,1],
                            [0,0,-a,0,0,0,a,0,-self.Γ31],
                            [0,L3E2,-b,0,0,0,0,a,0],
                            [-a,-b,L4E3,0,0,0,0,0,a],
                            [0,0,0,L5E4,0,-a,b,0,0],
                            [0,0,0,0,0,-b,0,b,-self.Γ32],
                            [0,0,0,-a,-b,L7E6,0,0,b],
                            [a,0,0,b,0,0,L7E8,0,-a],
                            [0,a,0,0,b,0,0,L7E9,-b],
                            [0,0,a,0,0,b,-a,-b,self.Γ3]]
                        ,dtype=complex)

        RHS = np.array([1,0,0,0,0,0,0,0,0,0]) # Tr(ρ)=1 and steady state solution
        
        # Solving the overdetermined system
        ρsol = np.dot(np.linalg.pinv(ρdot,1e-15),RHS) 
        
        # Hermitian density matrix condition
        ρ11 = np.real(ρsol[0])
        ρ22 = np.real(ρsol[4])
        ρ33 = np.real(ρsol[8])
        ρ12 = 0.5*(ρsol[1]+conj(ρsol[3]))
        ρ21 = 0.5*(ρsol[3]+conj(ρsol[1]))
        ρ13 = 0.5*(ρsol[2]+conj(ρsol[6]))
        ρ31 = 0.5*(ρsol[6]+conj(ρsol[2]))
        ρ23 = 0.5*(ρsol[5]+conj(ρsol[7]))
        ρ32 = 0.5*(ρsol[7]+conj(ρsol[5]))


        return [ρ11,ρ12,ρ13,ρ21,ρ22,ρ23,ρ31,ρ32,ρ33]



    def plotρ(self,savefig,showplot):

        """Plots the density matrix elements"""

        figsize = (12,5)

        s=r'$\rho_{11}$'
        
        fig, (ax1, ax2) = plt.subplots(1, 2,figsize=figsize)
        fig.suptitle(r'Density matrix element $\rho_{%s}$' % self.ρij[1:])
        ax1.set_xlabel(r'$\Delta_1$ (Rad Hz)')
        ax1.set_ylabel(r'Im[$\rho_{%s}$]' % self.ρij[1:])
        ax2.set_xlabel(r'$\Delta_1$ (Hz)')
        ax2.set_ylabel(r'Re[$\rho_{%s}$]' % self.ρij[1:])

        ax1.plot(self.Δrange,np.imag(self.ρΔid),'k-')
        ax2.plot(self.Δrange,np.real(self.ρΔid),'k-')

        ax1.grid()
        ax2.grid()
        fig.tight_layout()
        ax1.ticklabel_format(scilimits=(0,0))
        ax2.ticklabel_format(scilimits=(0,0))

        if savefig: plt.savefig("elements/" + self.ρij + ".png"); 
        if showplot: plt.show()


    def sweep_Δ1(self,ρid,fid):
        
        """Sweeps over a range of probe detunings and solves for the density elements"""

        self.ρij=ρijs[ρid]; self.ρid=ρid; self.fid=fid
        Δmax = (self.Ωp+self.Ωc) # This is an arbitrary choice
        self.Δrange=np.linspace(-Δmax,Δmax,self.fid+1)
        self.ρΔid = np.array([self.ρsol(Δ)[self.ρid] for Δ in self.Δrange])

    def __init__(self,Γ32,Γ31,γ3d,γ2d,Ωc,Ωp,Δ2):

        self.Γ32=Γ32; self.Γ31=Γ31; self.γ21=γ2d; self.Ωc=Ωc; self.Ωp=Ωp; self.Δ2=Δ2; self.γ3d=γ3d
        self.Γ3=Γ32+Γ31; self.γ31=self.Γ3+γ3d; self.γ32=self.γ31+self.γ21
