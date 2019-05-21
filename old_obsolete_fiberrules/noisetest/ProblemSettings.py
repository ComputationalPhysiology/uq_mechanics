from dolfin import *
import numpy as np
from MatStrainEnergyFunction2 import *
from fiberrules import *


#import cpp

parameters.allow_extrapolation = True
class problem():
	def __init__(self):

	  "initialization"
	  # Define mesh 
	  self.mesh = Mesh("refined_mesh.xml")#Mesh("rocio.xml")
	  
	  # Define normal vector 
          self.n = FacetNormal(self.mesh)

	  # Define body force
	  self.b = Constant((0.0, 0.0, 0.0))

          # Define function space
          self.V = VectorFunctionSpace(self.mesh, "CG", 1)

	  # Define fiber function space
	  self.fiber_space =   FunctionSpace(self.mesh,"CG",1) #self.V #FunctionSpace(self.mesh,"CG",1) 
          self.projection_space = FunctionSpace(self.mesh,"CG",1)


	  # Define trial and test functions
	  self.du = TrialFunction(self.V)              # Incremental displacement
	  self.u  = Function(self.V)                   # Displacement from previous iteration
	  self.v  = TestFunction(self.V)               # Test function
            
	  # Define kinematics
	  self.I = Identity(self.mesh.topology().dim())       # Identity tensor
          self.F = variable(self.I + grad(self.u))            # Deformation gradient

          # Create mesh functions over the cell facets
          self.subdomains = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1, self.mesh.domains() )


	  # Boundary markers, solutions and cases
	  self.markers = dict(base=10,rv=20,lv=30,epi=40,apex=50)


          # Define new measures associated with the exterior boundaries 
          self.ds = Measure("ds")[self.subdomains]

         # Define specific nodes as sub-domains 
          class Node1 (SubDomain):
                def inside(self,x,on_boundary):
                  return  near(x[0],0.0) and near(x[1], 1.274,0.9) and near(x[2], 1.892,0.9)



          class Node2 (SubDomain):
                def inside (self,x, on_boundary):
                  return  near(x[0],0.0) and near(x[1], -2.976,0.9) and near(x[2], 1.311,0.9)



          class Node3 (SubDomain):
                def inside(self,x,on_boundary):
                  return  near(x[0],0.0) and near(x[1], -1.844,0.9) and near(x[2], -2.191,0.9)



          class Node4 (SubDomain):
		def inside (self,x, on_boundary):
                  return  near(x[0],0.0) and near(x[1], 1.615,0.9) and near(x[2], -1.452,0.9)



	  class Node5 (SubDomain):
	        def inside (self,x,on_boundary):
		   return near(x[0],9.9,0.9) and near(x[1], 0.6056,0.9) and near(x[2], -0.1374,0.9)



          # 1- Basal plane: u_x = 0 for all nodes
          zero = Constant(0.0)
          bc0  = DirichletBC(self.V.sub(0), zero, self.subdomains, self.markers["base"],"topological")
	  self.bcs = bc0


          # 2- Basal plane: 4 nodes restrictions
          node1 = Node1()
          node2 = Node2()
          node3 = Node3()
          node4 = Node4()
  
	  #bc1  = DirichletBC(self.V.sub(2), zero, node1, method='pointwise')
          #bc2  = DirichletBC(self.V.sub(1), zero, node2, method='pointwise')
          #bc3  = DirichletBC(self.V.sub(2), zero, node3, method='pointwise')
          #bc4  = DirichletBC(self.V.sub(1), zero, node4, method='pointwise')
	  #self.bcs = [bc0,bc1,bc2,bc3, bc4]


          # 2- Apex: y-z restriction. Remove if Robin BC is running
	  node5 = Node5()
          bc5  = DirichletBC(self.V.sub(1), zero, node5, method='pointwise')
          bc6  = DirichletBC(self.V.sub(2), zero, node5, method='pointwise')
	  #self.bcs = [bc0,bc5,bc6]

        def update_params(self, b_ff = 6.6, b_xx = 4.0, b_fx = 2.6, CC = 1.1, K = 10.0, fiber_rotation_epi=50, fiber_rotation_endo=40,sheet_rotation_epi=65,sheet_rotation_endo=25,alpha_noise = 0.0, beta_noise = 0.0):
          self.b_ff = b_ff
          self.b_xx = b_xx
          self.b_fx = b_fx
          self.CC = CC
          self.K = K

          self.fiber_rotation_epi = fiber_rotation_epi
          self.fiber_rotaion_endo = fiber_rotation_endo
          self.sheet_rotation_epi = sheet_rotation_epi
          self.sheet_rotation_endo = sheet_rotation_endo

          self.alpha_noise = alpha_noise
          self.beta_noise = beta_noise
          

	  #  mesh and fiber function space
          mesh = self.mesh
          VV = self.fiber_space

          #print  b_ff , b_xx , b_fx, CC , K , fiber_rotation_epi, fiber_rotation_endo,sheet_rotation_epi,sheet_rotation_endo,alpha_noise , beta_noise 

          # Define fibre orientation field
          f0, s0, n0 = dolfin_fiberrules(mesh,VV, fiber_rotation_epi, fiber_rotation_endo,sheet_rotation_epi,sheet_rotation_endo, alpha_noise, beta_noise)

          self.f0 = f0#Constant([1.0,0.0,0.0])
          self.s0 = s0#Constant([0.0, 1.0, 0.0])
          self.n0 = n0#Constant([0.0, 0.0, 1.0])

          return





        def update_randomfield(self,  b_ff = 6.6, b_xx = 4.0, b_fx = 2.6, CC = 1.1, K = 10.0,  fiber_rotation_epi=50, fiber_rotation_endo=40,sheet_rotation_epi=65,sheet_rotation_endo=25, alpha_noise = 0.0, beta_noise = 0.0):

          self.b_ff = b_ff
          self.b_xx = b_xx
          self.b_fx = b_fx
          self.CC = CC
          self.K = K

          self.fiber_rotation_epi = fiber_rotation_epi
          self.fiber_rotaion_endo = fiber_rotation_endo
          self.sheet_rotation_epi = sheet_rotation_epi
          self.sheet_rotation_endo = sheet_rotation_endo
          self.alpha_noise = alpha_noise
          self.beta_noise = beta_noise

          #  mesh and fiber function space
          mesh = self.mesh
          VV = self.fiber_space
          print  b_ff , b_xx , b_fx, CC , K , fiber_rotation_epi, fiber_rotation_endo,sheet_rotation_epi,sheet_rotation_endo, beta_noise


          # Define Fibre orientation Field
          f0,s0,n0 =  dolfin_fiberrules(mesh,VV, fiber_rotation_epi, fiber_rotation_endo,sheet_rotation_epi,sheet_rotation_endo,alpha_noise,beta_noise)
          self.f0 = f0#Constant([1.0,0.0,0.0])
          self.s0 = s0#Constant([0.0, 1.0, 0.0])
          self.n0 = n0#Constant([0.0, 0.0, 1.0])




	def run (self,pressure):
      
	  mesh = self.mesh #mesh
          n = self.n       #normal
	  ds = self.ds     #subdomains 
	  markers = self.markers #subdomain regions  
	  b = self.b       #body force
    	  V = self.V       #function space
	  bcs = self.bcs   #boundary conditions

          #trial and test functions
	  du = self.du               # Incremental displacement
	  u  = self.u                # Displacement from previous iteration
	  v  = self.v                # Test function
          projection_space = self.projection_space           

	  # Define kinematics
	  I = self.I                 # Identity tensor
          F = self.F                 # Deformation gradient


          # parameters
	  b_ff = self.b_ff
	  b_xx = self.b_xx
	  b_fx = self.b_fx
          CC = self.CC
	  K = self.K

	  f0 = self.f0
          n0 = self.n0
          s0 = self.s0

      	  # Define strain energy density
	  material =  StrainEnergyFunction(I,F) #StrainEnergyFunction(F)
          Psi_iso = material.psi_iso(b_ff, b_xx, b_fx, CC, K, f0, n0, s0)
          Psi_vol = material.psi_vol(K)
	  Psii = Psi_iso + Psi_vol

	  # Define the elastic response of the material 
	  # First Piola-Kircchoff tensor
	  P = diff(Psii, F)
	  # Second Piola-Kircchoff stress tensor
          #S = 2.0*diff(Psi, C)

	  # Define  variational formulation first term
	  Psi = inner(P, grad(v))*dx


          #Pressure
    	  p = Constant(-pressure)
      		
	  # Robin Boundary Condition term
          k = Constant(1.0) #kPa
	  rbc  = k*inner(u,v)*ds(markers["base"])

    	  # Total Potential Energy 
    	  Pi = Psi  + rbc - inner(p*n, v)*ds(markers["lv"]) - inner(b,v)*dx 
        
          # Differenciate to get the Jacobian of Pi 
    	  Jac  = derivative(Pi, u, du)

    	  # Solve the boundary condition problem
    	  solve (Pi == 0, u , bcs, J=Jac, solver_parameters={"newton_solver":{"error_on_nonconvergence": False, "maximum_iterations": 50}})
        
          #plot(u, interactive = True)

	  #quantities of interest
	  #increase in inner Cavity Volume
	  dV = self.icavity_volume(u)

          #lengthening of the Ventricle
    	  dL = self.lenghtening(u)
 
    	  #increase in inner radius
    	  dR_endo = self.inner_radius(u)

          #increase in ext. radius
	  dR_epi = self.ext_radius(u)

  	  # Calculate decrease thickness ventricle wall
    	  dT = self.ventwall_thickness(dR_epi,dR_endo)

    	  # Calculate  twisting at base
	  dTwist = self.twisting_angle(u)

          # Calculate projecton of the Second Piola-Kirchhoff stress tensor onto the fiber axis direction
          S = material.SPK_stresstensor(b_ff, b_xx, b_fx, CC, K, f0, n0, s0,u)
          #TV = TensorFunctionSpace(mesh, 'DG',0)
          #Sp = project(S,TV)
          #plot(Sp[0,0],interactive=True)
          Sff = material.SPK_projection(S,f0)
          Sffp = project(Sff,projection_space)
          #plot(Sffp,interactive=True)
          #print Sffp.vector().array().shape
          #print u.vector().array().shape

          # Output list
 	  output = [dV,dL, dR_endo, dR_epi, dT, dTwist,  u.vector().array(), Sffp.vector().array()]
          return output #u.vector().array()


	def icavity_volume(self, u):
	  # Calculate increase in Cavity Volume
	  n = self.n
          ds = self.ds
          mesh= self.mesh
          subdomains = self.subdomains
          dV =  assemble(inner(u,n) * ds(self.markers["lv"] , subdomain_data = self.subdomains))
    	  return dV	

	def lenghtening(self, u):
           # Calculate lengthening of the Ventricle
           return u(9.90194,0.695668,-0.137405)[0]

	def inner_radius(self,u):
	  # Calculate increase in inner radius
	  ir2 = u(0.0, 1.274,1.892)[1]**2 + u(0.0, 1.274,1.892)[2]**2
          return sqrt(ir2)

	def ext_radius(self,u,):
          # Calculate increase in ext. radius
	  er2 = u(0.0, 1.867,2.6467)[1]**2 + u(0.0, 1.867,2.6467)[2]**2
          return sqrt(er2)

	def ventwall_thickness(self, dr_epi,dr_endo):
  	  # Calculate decrease thickness ventricle wall
	  return dr_epi - dr_endo	  

	def twisting_angle(self,u):
	  # Calculate  twisting at base
	  dy = u(0.0,1.867,2.6467)[1]
    	  dz = u(0.0,1.867,2.6467)[2]
	  dtheta = np.arctan2(dy,dz)
	  return dtheta

