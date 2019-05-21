from dolfin import *
import numpy as np
from MatStrainEnergyFunction2 import *
from fiberrules import *
from meshutils import *

parameters.allow_extrapolation = True
class problem():
    def __init__(self):
        # Define mesh 
        self.mesh = Mesh("refined_mesh.xml")#Mesh("rocio.xml")
        # Define normal vector 
        self.n = FacetNormal(self.mesh)

        # Define body force
        self.b = Constant((0.0, 0.0, 0.0))

        # Define function space
        degree=2
        self.V = VectorFunctionSpace(self.mesh, "CG", degree)

        # Define trial and test functions
        self.du = TrialFunction(self.V)              # Incremental displacement
        self.u  = Function(self.V)                   # Displacement from previous iteration
        self.v  = TestFunction(self.V)               # Test function
            
        # Define kinematics
        self.I = Identity(self.mesh.topology().dim())       # Identity tensor
        self.F = variable(self.I + grad(self.u))            # Deformation gradient

        # Create mesh functions over the cell facets
        self.subdomains = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1, self.mesh.domains() )
        
        # Create mesh functions over the domains
        self.domain = self.mesh.ufl_domain()

        # Define fiber function space
        self.fiber_space =   FunctionSpace(self.mesh,"CG",degree-1) #self.V #FunctionSpace(self.mesh,"CG",2) 
        self.projection_space = FunctionSpace(self.mesh,"CG",degree)

        # Boundary markers, solutions and cases
        self.markers = dict(base=10,rv=20,lv=30,epi=40,apex=50)


        # Define new measures associated with the exterior boundaries 
        #self.ds = Measure("ds")[self.subdomains]
        self.ds = ds(subdomain_data = self.subdomains)


        # 1- Basal plane: u_x = 0 for all nodes
        zero = Constant(0.0)
        bc0  = DirichletBC(self.V.sub(0), zero, self.subdomains, self.markers["base"],"topological")
        self.bcs = bc0



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
        #f0, s0, n0 = dolfin_fiberrules(mesh,VV, fiber_rotation_epi, fiber_rotation_endo,sheet_rotation_epi,sheet_rotation_endo, alpha_noise, beta_noise)
        f0, s0, n0 = dolfin_fiberrules(mesh,VV, fiber_rotation_epi, fiber_rotation_endo,sheet_rotation_epi,sheet_rotation_endo)

        
        self.f0 = f0#Constant([1.0,0.0,0.0])
        self.s0 = s0#Constant([0.0, 1.0, 0.0])
        self.n0 = n0#Constant([0.0, 0.0, 1.0])






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
        print(alpha_noise)
        print(beta_noise)

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
        solve (Pi == 0, u , bcs, J=Jac, solver_parameters={"newton_solver":{'linear_solver':'mumps',"error_on_nonconvergence": False, "maximum_iterations": 15}})
        

        #plot(u, interactive = True)

        if np.isnan([u(0,0,0)[0]]):
            output = np.empty(8)
            output[:] = np.NAN
        else:

            #quantities of interest
            #increase in inner Cavity Volume
            #dV = self.icavity_volume(u)
            # Total inner volume
            dV = self.inner_volume(u)

            #lengthening of the Ventricle
            #dL = self.lenghtening(u)
            # total le
            dL_endo = 0#self.ventricle_axial_length("lv",u)
            dL_epi = 0#self.ventricle_axial_length("epi",u)
            dL = 0#self.ventwall_thickness(dL_epi,dL_endo)


    	    #increase in inner radius
    	    #dR_endo = self.inner_radius(u)
            dR_endo= 0#self.average_radius("lv",u)
     
            #increase in ext. radius
            #dR_epi = self.ext_radius(u)
            dR_epi = 0#self.average_radius("epi",u)

            # Calculate decrease thickness ventricle wall
            dT = 0# self.ventwall_thickness(dR_epi,dR_endo)

            # Calculate rotation at apex
            dTwist = 0#self.apex_rotation(u)

    	    # Calculate  twisting at base
            #dTwist = self.twisting_angle(u)

            # Calculate the wall volume
            dWV = 0#self.wall_volume(u)

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
            output = [dV,dL_endo,dL_epi,dL, dR_endo, dR_epi, dT, dTwist, dWV , u.vector().array(), Sffp.vector().array()]

            return output #u.vector().array()


    def inner_volume(self, u = Constant((0.0, 0.0, 0.0))) :
        """
        Compute the inner volume of the cavity for a given displacement u.
        """

        # current inner volume
        domain = self.domain
        X = SpatialCoordinate(self.mesh)
        N = FacetNormal(self.mesh)

        x = X + u
        F = grad(x)
        n = cofac(F) * N

        ds_endo = ds(self.markers["lv"], subdomain_data = self.subdomains)
        
        Vendo_form = -1/float(3) * inner(x, n) * ds_endo

        return assemble(Vendo_form)




    def ventricle_axial_length(self, surface, u = None) :
        """
        Compute the length of the ventricle by intersecting it with the
        line directed as its axis. 
        surface = "lv" or "epi"
        """
        # restrict on boundary
        marker = self.markers[surface]
        meshb = boundary_mesh(self.mesh, self.subdomains, marker)
        if u :
            ub = restrict_on_boundary(meshb, u)
        else :
            ub = None

        # vtk grid
        grid = dolfin2vtk(meshb, ub)
        lingeo = linear_grid_filter(grid, 1)

        long_axis = 0

        return lingeo.GetOutput().GetBounds()[2*long_axis + 1]




    def average_radius(self, surface, u = None) :
        """
        Compute average radius of an axial section of the ventricle.
        """


        # restrict on boundary
        marker = self.markers[surface]
        meshb = boundary_mesh(self.mesh, self.subdomains, marker)
        if u :
            ub = restrict_on_boundary(meshb, u)
        else :
            ub = None

        # vtk grid
        grid = dolfin2vtk(meshb, ub)
        lingeo = linear_grid_filter(grid, 1)

        # cut plane
        long_axis = 0
        origin = [ 0.0, 0.0, 0.0 ]
        normal = [ 0.0, 0.0, 0.0 ]
        origin[long_axis] = 1.0 #cut_quote
        normal[long_axis] = 1.0
        plane = vtkPlane()
        plane.SetOrigin(*origin)
        plane.SetNormal(*origin)

        plane_cut = vtkCutter()
        if vtk6:
            plane_cut.SetInputData(lingeo.GetOutput())
        else:
            plane_cut.SetInput(lingeo.GetOutput())
            
        plane_cut.SetCutFunction(plane)
        plane_cut.Update()

        rsum = 0.0
        pts = plane_cut.GetOutput().GetPoints()
        for i in xrange(0, pts.GetNumberOfPoints()) :
            rsum += np.linalg.norm(pts.GetPoint(i)[1:])

        return rsum / pts.GetNumberOfPoints()


    def wall_volume(self, u = Constant((0.0, 0.0, 0.0))) :
        """
        Compute the volume of the wall.
        """

        X = SpatialCoordinate(self.mesh)
        x = X - u #original was X + u
        F = grad(x)
        J = det(F)

        V_form = J * dx(domain = self.mesh)

        return assemble(V_form)

    def apex_rotation(self, u = Constant((0.0, 0.0, 0.0))):
        """ Compute the rotation angle of the apex
        Apex Coordinates: 9.89854832122, 0.284472435693, -0.224218173585 """
        x0 = 9.898548
        y0 = 0.284472
        z0 = -0.22422

        # Apex coordinates after displacement
        x = x0 + u(x0,y0,z0)[0]
        y = y0 + u(x0,y0,z0)[1]
        z = z0 + u(x0,y0,z0)[2]

        sprod =  y0*y + z0*z
        vec0_norm = sqrt(y0*y0+z0*z0)
        vec_norm =  sqrt(y*y + z*z)

        costheta = sprod/(vec0_norm*vec_norm)
        
        theta = np.arccos(costheta)
        return theta




    def icavity_volume(self, u):
        # Calculate increase in Cavity Volume
        n = self.n
        mesh = self.mesh
        ds = self.ds
        subdomains = self.subdomains
        dV =  assemble(inner(u,n) * ds(self.markers["lv"] , subdomain_data = self.subdomains))
        return dV


    def ventwall_thickness(self, dr_epi,dr_endo):
        # Calculate decrease thickness ventricle wall
	    return dr_epi - dr_endo	  

