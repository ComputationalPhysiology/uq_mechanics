from dolfin import *



class StrainEnergyFunction():
  # Define the material: 

  def __init__(self, I,F):
     self.I = I
     self.F = F
         
  def psi_iso(self, b_ff, b_xx, b_fx, CC, K, f0, n0, s0):

    F = self.F
    
    # Right Cauchy-Green tensor
    C = F.T*F
    # Green-Lagrange strain tensor
    I = self.I#Identity(mesh.topology().dim())
    E = 0.5 *  (C - I)
    

    E_ff = dot(f0, E*f0) #WARNING: inner or dot??!!!
    E_nn = dot(n0, E*n0)
    E_ss = dot(s0, E*s0)
    E_sn = dot(s0, E*n0)
    E_ns = dot(n0, E*s0)
    E_fn = dot(f0, E*n0)
    E_nf = dot(n0, E*f0)
    E_fs = dot(f0, E*s0)
    E_sf = dot(s0, E*f0)


    W1 = b_ff * E_ff**2
    W2 = b_xx * ( E_nn**2 + E_ss**2 + E_sn**2 + E_ns**2)
    W3 = b_fx * ( E_fn**2 + E_nf**2  + E_fs**2 + E_sf**2)
    W = W1 + W2 + W3


    isochoric = 0.5 * CC *(exp(W) - 1.0)
    return isochoric




  def psi_iso_invariants(self, b_ff, b_xx, b_fx, CC, K, f0, n0, s0):
    # 1 - Isochoric part

     F = self.F
     
     # Right Cauchy-Green tensor
     C = F.T*F

     # Isotropic Invariants
     I1 = tr(C)
     I2 = 0.5*(tr(C)**2 - tr(C*C))

     # Redefine constants
     a_1 = 0.25 * (b_ff - b_xx)
     a_2 = 0.25 * b_xx
     a_3 = 0.5 * (b_fx-b_xx)

     C_ff = dot(f0, C*f0) #WARNING: inner or dot??!!!
     C_fn = dot(f0, C*n0)
     C_fs = dot(f0, C*s0)


     W1 = a_1 * ( C_ff - 1.)**2
     W2 = a_2 * ( I1**2 - 2.*(I1+I2) + 3.)
     W3 = a_3 * (C_fn**2 + C_fs**2)
     W = W1 + W2 + W3

     #isochoric strain energy term
     isochoric = 0.5 * CC *(exp(W) - 1.0)

     return isochoric


  def psi_vol(self, K):
    # 2 - Volumetric part (almost incompressible material)
     F  = self.F

     J = det (F)
     volumetric = K* (J * ln(J) - J + 1. )

     return volumetric

  def SPK_stresstensor(self, b_ff, b_xx, b_fx, CC, K, f0, n0, s0, u):
     
    I = self.I
    F = I + grad(u)
    J = det(F)

    # Right Cauchy-Green tensor
    C = F.T*F
    # Green-Lagrange strain tensor
    E = 0.5 *  (C - I)


    E_ff = dot(f0, E*f0) #WARNING: inner or dot??!!!
    E_nn = dot(n0, E*n0)
    E_ss = dot(s0, E*s0)
    E_sn = dot(s0, E*n0)
    E_ns = dot(n0, E*s0)
    E_fn = dot(f0, E*n0)
    E_nf = dot(n0, E*f0)
    E_fs = dot(f0, E*s0)
    E_sf = dot(s0, E*f0)

    # Exponent of Strain Energy function Psi
    W1 = b_ff * E_ff**2
    W2 = b_xx * ( E_nn**2 + E_ss**2 + E_sn**2 + E_ns**2)
    W3 = b_fx * ( E_fn**2 + E_nf**2  + E_fs**2 + E_sf**2)
    W = W1 + W2 + W3

    # Derivative of Psi_iso respect to C
    dEff_dC = outer(f0,f0) #f0*f0.T
    dEss_dC = outer(s0,s0) #s0*s0.T
    dEnn_dC = outer(n0,n0) #n0*n0.T
    dEsn_dC = outer(s0,n0) #s0*n0.T
    dEns_dC = outer(n0,s0) #n0*s0.T
    dEfn_dC = outer(f0,n0) #f0*n0.T
    dEnf_dC = outer(n0,f0) #n0*f0.T
    dEfs_dC = outer(f0,s0) #f0*s0.T
    dEsf_dC = outer(s0,f0) #s0*f0.T

    dWdC1 = 2.0*b_ff*(E_ff*dEff_dC)
    dWdC2 = 2.0*b_xx*(E_ss*dEss_dC + E_nn*dEnn_dC+  E_sn*dEsn_dC + E_ns*dEns_dC) 
    dWdC3 = 2.0*b_fx*(E_fn*dEfn_dC + E_nf*dEnf_dC+  E_fs*dEfs_dC + E_sf*dEsf_dC)
    dWdC = dWdC1 + dWdC2 + dWdC3

    dPsidC_iso = 0.5*CC*dWdC*exp(W)

    # Derivative of Psi_vol respect to C
    dPsidJ = K * ln(J) 
    dJdC = J*inv(C)
    dPsidC_vol = dPsidJ * dJdC 

    # Second Piola-Kirchhoff stress tensor
    S = 2.0*(dPsidC_iso +  dPsidC_vol)
    
    return S


  def SPK_projection(self,S,a0):
     return  inner(a0,S*a0)


