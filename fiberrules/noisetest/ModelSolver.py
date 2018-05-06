from dolfin import *
import numpy as np
from inputs2 import *
from FibreField2 import *
from MatStrainEnergyFunction2 import *

class problem():
	def __init__(self):

	  "initialization"
          # Import geometry parameters
	  self.geoparam = inputs.geometry_parameters()
	  a_endo = self.geoparam['a_endo']
	  b_endo = self.geoparam['b_endo']
	  c_endo = self.geoparam['c_endo']
	  a_epi = self.geoparam['a_epi']
	  b_epi = self.geoparam['b_epi']
	  c_epi = self.geoparam['c_epi']


	  # Define mesh 
	  self.mesh = Mesh("leftvent_art.xml")
	  
	  # Define normal vector 
          self.n = FacetNormal(self.mesh)

	  # Define body force
	  self.b = Constant((0.0, 0.0, 0.0))

          # Define function spaces
          self.V = VectorFunctionSpace(self.mesh, "CG", 1)

	  # Define trial and test functions
	  self.du = TrialFunction(self.V)              # Incremental displacement
	  self.u  = Function(self.V)                   # Displacement from previous iteration
	  self.v  = TestFunction(self.V)               # Test function

	  # Define kinematics
	  self.I = Identity(self.mesh.topology().dim())  # Identity tensor
          self.F = variable(self.I + grad(self.u))            # Deformation gradient
	  
	  # Define exterior boundary subdomains
	  class Basal_Plane(SubDomain):
          	def inside(self,x, on_boundary):
                  return on_boundary and near(x[0], 0.0)

	  class Endocardium(SubDomain):
          	def inside(self, x, on_boundary):
                  r_endo = x[1]*x[1]/a_endo**2 + x[2]*x[2]/b_endo**2 + x[0]*x[0]/c_endo**2
                  return on_boundary  and near(r_endo, 1.0, 0.09)

	  class Epicardium(SubDomain):
          	def inside(self,x, on_boundary):
                  r_epi = x[1]*x[1]/a_epi**2 + x[2]*x[2]/b_epi**2 + x[0]*x[0]/c_epi**2
                  return on_boundary  and near(r_epi, 1.0, 0.01)


	  # Initialize sub-domain instances
	  base = Basal_Plane()
	  endo = Endocardium()
	  epi =  Epicardium()

          # Create mesh functions over the cell facets
          self.subdomains = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1 )

          # Mark all facets as subdomain 3
          self.subdomains.set_all(3)
          # Mark Base boundary as subdomain 0
          base.mark(self.subdomains, 0)
          # Mark endocardium boundary as subdomain 1
          endo.mark(self.subdomains, 1)
          # Mark Epicardium boundary as subdomain 2
          epi.mark(self.subdomains, 2)


          # Define new measures associated with the exterior boundaries 
          self.ds = Measure("ds")[self.subdomains]
	  
	  file = File("subdomains.pvd")
	  file << self.subdomains

	  # Define specific nodes as sub-domains 
	  class Node1 (SubDomain):
          	def inside(self,x,on_boundary):
          	  return  near(x[0],0.0) and near(x[1], a_endo) and near(x[2], 0.0)



	  class Node2 (SubDomain):
          	def inside (self,x, on_boundary):
                  return  near(x[0],0.0) and near(x[1], 0.0) and near(x[2], b_endo)



	  class Node3 (SubDomain):
          	def inside(self,x,on_boundary):
          	  return  near(x[0],0.0) and near(x[1], -a_endo) and near(x[2], 0.0)



	  class Node4 (SubDomain):
          	def inside (self,x, on_boundary):
                  return  near(x[0],0.0) and near(x[1], 0.0) and near(x[2], -b_endo)



          # Define boundary conditions 

          # 1- Basal plane: u_x = 0 for all nodes
          zero = Constant(0.0)
          bc0  = DirichletBC(self.V.sub(0), zero, base)

          # 2- Basal plane: 4 nodes restrictions
	  node1 = Node1()
	  node2 = Node2()
          node3 = Node3()
	  node4 = Node4()

	  bc1  = DirichletBC(self.V.sub(2), zero, node1, method='pointwise')
          bc2  = DirichletBC(self.V.sub(1), zero, node2, method='pointwise')
          bc3  = DirichletBC(self.V.sub(2), zero, node3, method='pointwise')
          bc4  = DirichletBC(self.V.sub(1), zero, node4, method='pointwise')

          self.bcs = [bc0,bc1,bc2,bc3, bc4]

	  



	  # random parameters default values
	  #self.alpha =  60.0
	  #self.b_ff = 1#b_ff
	  #self.b_xx = 1#b_xx
	  #self.b_fx = 1##b_fx
          #self.CC = 1#CC
	  #self.K = 1#K



	def update_params(self,alpha,b_ff, b_xx, b_fx, CC, K ):
	  self.alpha =  alpha
	  self.b_ff = b_ff
	  self.b_xx = b_xx
	  self.b_fx = b_fx
          self.CC = CC
	  self.K = K
	  
	  # Define fibre orientation field
	  fiber =  FibreField(self.mesh, self.geoparam, 60,-60)

	  self.f0 = fiber.f0#Constant([1.0,0.0,0.0])
	  self.s0 = fiber.s0#Constant([0.0, 1.0, 0.0])
          self.n0 = fiber.n0#Constant([0.0, 0.0, 1.0])

	  return


	def run (self,pressure):

	  self.geoparam = inputs.geometry_parameters()
	  a_endo = self.geoparam['a_endo']
	  b_endo = self.geoparam['b_endo']
	  c_endo = self.geoparam['c_endo']
	  a_epi = self.geoparam['a_epi']
	  b_epi = self.geoparam['b_epi']
	  c_epi = self.geoparam['c_epi']
       
	  mesh = self.mesh #mesh
          n = self.n       #normal
	  ds = self.ds     #subdomains   
	  b = self.b       #body force
    	  V = self.V       #function space
	  bcs = self.bcs   #boundary conditions

          #trial and test functions
	  du = self.du               # Incremental displacement
	  u  = self.u                # Displacement from previous iteration
	  v  = self.v                # Test function

	  #Define kinematics
	  I = self.I                 # Identity tensor
          F = self.F                 # Deformation gradient

          #random parameters
	  alpha =  self.alpha
	  b_ff = self.b_ff
	  b_xx = self.b_xx
	  b_fx = self.b_fx
          CC = self.CC
	  K = self.K

	  f0 = self.f0
          n0 = self.n0
          s0 = self.s0

      	  # Define strain energy density
	  material = StrainEnergyFunction(F)
          Psi_iso = material.psi_iso(b_ff, b_xx, b_fx, CC, K, f0, n0, s0)
          Psi_vol = material.psi_vol(K)
	  Psi = Psi_iso + Psi_vol

	  # Define the elastic response of the material 
	  # First Piola-Kircchoff tensor
	  P = diff(Psi, F)

	  # Define  variational formulation first term
	  Psi = inner(P,grad(v))*dx

    	  p = Constant(-pressure)

      		
    	  # Total Potential Energy 
    	  Pi = Psi  - inner(p*n, v)*ds(1) - inner(b,v)*dx
        
          # Differenciate to get the Jacobian of Pi 
    	  Jac  = derivative(Pi, u, du)

    	  # Solve the boundary condition problem
    	  solve (Pi == 0, u , bcs, J=Jac )
        
          plot(u, interactive = True)

	  #quantities of interest
	  #increase in inner Cavity Volume
	  dV = self.icavity_volume(u)

          #lengthening of the Ventricle
    	  dL = self.lenghtening(u, c_epi)
 
    	  #increase in inner radius
    	  dR_endo = self.inner_radius(u,b_endo)

          #increase in ext. radius
	  dR_epi = self.ext_radius(u,b_epi)

  	  # Calculate decrease thickness ventricle wall
    	  dT = self.ventwall_thickness(dR_epi,dR_endo)

    	  # Calculate  twisting at base
	  dTwist = self.twisting_angle(u, b_epi)

	  output = [dV,dL, dR_endo, dR_epi, dT, dTwist]
          return u#output#u.vector().array()


	def icavity_volume(self, u):
	  # Calculate increase in Cavity Volume
	  n = self.n
          ds = self.ds
          mesh= self.mesh
          subdomains =self.subdomains
          dV  = assemble(dot(u,n) * ds(1, domain = mesh, subdomain_data = subdomains))
    	  return dV	

	def lenghtening(self, u, coord):
           # Calculate lengthening of the Ventricle
           return u(coord,0.0,0.0)[0]

	def inner_radius(self,u, coord):
	  # Calculate increase in inner radius
          return u(0.0,0.0,coord)[2]

	def ext_radius(self,u, coord):
          # Calculate increase in ext. radius
          return u(0.0,0.0,coord)[2]

	def ventwall_thickness(self, dr_epi,dr_endo):
  	  # Calculate decrease thickness ventricle wall
	  return dr_epi - dr_endo	  

	def twisting_angle(self,u, coord):
	  # Calculate  twisting at base
	  dy = u(0.0,0.0,coord)[1]
    	  dz = u(0.0,0.0,coord)[2]
	  dtheta = np.arctan2(dy,dz)
	  return dtheta

