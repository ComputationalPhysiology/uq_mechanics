__all__ = ["dolfin_fiberrules", "dolfin_to_carpfile", "dolfin_to_vtk", \
           "scalar_laplacians"]
import dolfin as d
import numpy as np
from kle import KLE
from covariance import CovarianceMatrix
from kernel import Matern
import fiberrules as cpp


def dolfin_fiberrules(mesh, fiber_space=None,
                      fiber_rotation_epi=50,  # 50 
                      fiber_rotation_endo=40, # 40
                      sheet_rotation_epi=65,  # 65
                      sheet_rotation_endo=25, # 25
                      alpha_noise=0.0 ,beta_noise=0.0):
    """
    Create fiber, cross fibers and sheet directions

    Arguments
    ---------
    mesh : dolfin.Mesh
       A dolfin mesh with marked boundaries:
       base = 10, rv = 20, lv = 30, epi = 40
       The base is assumed placed at x=0
    fiber_space : dolfin.FunctionSpace (optional)
       Determines for what space the fibers should be calculated for.
    fiber_rotation_epi : float (optional)
       Fiber rotation angle on the endocardial surfaces.
    fiber_rotation_endo : float (optional)
       Fiber rotation angle on the epicardial surfaces.
    sheet_rotation_epi : float (optional)
       Sheet rotation angle on the endocardial surfaces.
    sheet_rotation_endo : float (optional)
       Sheet rotation angle on the epicardial surfaces.
    """
    #import cpp
    
    if not isinstance(mesh, d.Mesh):
        raise TypeError("Expected a dolfin.Mesh as the mesh argument.")

    # Default fiber space is P1
    fiber_space = fiber_space or d.FunctionSpace(mesh, "P", 1)

    # Create scalar laplacian solutions
    d.info("Calculating scalar fields")
    scalar_solutions = scalar_laplacians(mesh)

    # Create gradients
    d.info("\nCalculating gradients")
    data = project_gradients(mesh, fiber_space, scalar_solutions)
 
    # Assign the fiber and sheet rotations
    data.fiber_rotation_epi = fiber_rotation_epi
    data.fiber_rotation_endo = fiber_rotation_endo
    data.sheet_rotation_epi = sheet_rotation_epi
    data.sheet_rotation_endo = sheet_rotation_endo



    #Gaussian field with correlation length l
    #p , l = np.inf, 2
    #number of terms in expansion
    #k = 10




    #kernel = Matern(p = p,l = l)    #0,1Exponential covariance kernel
    #kle  = KLE(mesh, kernel, verbose = True)
    #kle.compute_eigendecomposition(k = k)#20 by default

    #Generate realizations
    #noise = kle.realizations()
    #x = np.zeros(len(noise))
    #noise = np.array(x)
    #noise[0:len(noise)/2] = 30
    ##print y
    #print "noise array", noise



    #V = d.FunctionSpace(mesh, "CG", 1)
    #random_field =d.Function(V)
    #random_f = np.zeros((mesh.num_vertices()),dtype =float) 
    #random_f = noise
    #random_field.vector()[:] = random_f[d.dof_to_vertex_map(V)]
    #d.plot(random_field, mesh, interactive =True)

    ##random_f = noise
    
    #x = np.zeros(len(noise))
    #y = np.array(x)
    #y[0:len(noise)/2] = 30
    ##print y
    #print "noise array", y 

    # Check noise
    if np.isscalar(alpha_noise):
       alpha_noise = np.zeros(mesh.num_vertices(),dtype =float)

    if np.isscalar(beta_noise):
       beta_noise = np.zeros(mesh.num_vertices(),dtype =float)


    # Call the fiber sheet generation
    cpp.computeFiberSheetSystem(data,alpha_noise[d.dof_to_vertex_map(fiber_space)],beta_noise[d.dof_to_vertex_map(fiber_space)])
    #cpp.computeFiberSheetSystem(data,noise)



    # Create output Functions
    Vv = d.VectorFunctionSpace(mesh, fiber_space.ufl_element().family(), \
                               fiber_space.ufl_element().degree())
    V = fiber_space

    fiber_sheet_tensor = data.fiber_sheet_tensor
    fiber_components = []
    scalar_size = fiber_sheet_tensor.size/9
    indices = np.zeros(scalar_size*3, dtype="L")
    indices[0::3] = np.arange(scalar_size, dtype="L")*9   # x
    indices[1::3] = np.arange(scalar_size, dtype="L")*9+3 # y
    indices[2::3] = np.arange(scalar_size, dtype="L")*9+6 # z
    
    for ind, name in enumerate(["f0", "n0", "s0"]):
        component = d.Function(Vv, name=name)
        
        # Sort the fibers and sheets in dolfin degrees of freedom (dofs)
        component.vector()[:] = fiber_sheet_tensor[indices]
        fiber_components.append(component)
        indices += 1

    # Make sheet the last component
    fiber_components = [fiber_components[0]]+fiber_components[-1:0:-1]

    return fiber_components

def project_gradients(mesh, fiber_space, scalar_solutions):
    """
    Calculate the gradients using projections
    """
    #import cpp
    
    scalar_arrays = {}
    gradient_arrays = {}
    gradients = {}

    # If linear lagrange for fibers we do not need to interpolate the scalar solutions
    if fiber_space.ufl_element().family()=="Lagrange" and \
       fiber_space.ufl_element().degree()==1:
        lagrange = True
        scalar_interp = scalar_solutions["epi"]
    else:
        # Assume 3D and create a VectorSpace of the fiber space
        scalar_interp = d.Function(fiber_space)
        lagrange = False

    Vv = d.VectorFunctionSpace(mesh, fiber_space.ufl_element().family(), \
                               fiber_space.ufl_element().degree())
    # Pick one local solution
    local_size = scalar_interp.vector().local_size()

    # Generate the FiberRulesData object which will be passed to FiberSheet generation
    data = cpp.FiberRulesData(local_size)
    for case, scalar_solution in scalar_solutions.items():
        gradient = d.project(d.grad(scalar_solution), Vv)

        if not lagrange:
            scalar_solution = d.interpolate(scalar_solution, fiber_space)

        # Add scalar data
        if case != "apex":
            data[case+"_scalar"] = scalar_solution.vector().array()

        # Add gradient data
        data[case+"_gradient"] = gradient.vector().array()
            
    # Return data
    return data

def scalar_laplacians(mesh):
    """
    Calculate the laplacians needed by fiberrule algorithms
    
    Arguments
    ---------
    mesh : dolfin.Mesh
       A dolfin mesh with marked boundaries:
       base = 10, rv = 20, lv = 30, epi = 40
       The base is assumed placed at x=0
    
    """
    
    if not isinstance(mesh, d.Mesh):
        raise TypeError("Expected a dolfin.Mesh as the mesh argument.")

    # Init connectivities
    mesh.init(2)
    facet_markers = d.MeshFunction("size_t", mesh, 2, mesh.domains())

    # Boundary markers, solutions and cases
    markers = dict(base=10, rv=20, lv=30, epi=40, apex=50)

    # Solver parameters
    solver_param=dict(solver_parameters=dict(
        preconditioner="ml_amg" if d.has_krylov_solver_preconditioner("ml_amg") \
        else "default", linear_solver="gmres"))
    
    cases = ["rv", "lv", "epi"]
    boundaries = cases + ["base"]

    # Check that all boundary faces are marked
    num_boundary_facets = d.BoundaryMesh(mesh, "exterior").num_cells()

    if num_boundary_facets != sum(np.sum(\
        facet_markers.array()==markers[boundary])\
                                  for boundary in boundaries):
        d.error("Not all boundary faces are marked correctly. Make sure all "\
                "boundary facets are marked as: base = 10, rv = 20, lv = 30, "\
                "epi = 40.")
    
    # Coords and cells
    coords = mesh.coordinates()
    cells_info = mesh.cells()
    
    # Find apex by solving a laplacian with base solution = 0
    # Create Base variational problem
    V = d.FunctionSpace(mesh, "CG", 1)
    
    u = d.TrialFunction(V)
    v = d.TestFunction(V)
    
    a = d.dot(d.grad(u), d.grad(v))*d.dx
    L = v*d.Constant(1)*d.dx
    
    DBC_10 = d.DirichletBC(V, 1, facet_markers, markers["base"], "topological")

    # Create solutions
    solutions = dict((what, d.Function(V)) for what in markers if what != "base")

    d.solve(a==L, solutions["apex"], DBC_10,
            solver_parameters={"linear_solver": "gmres"})

    apex_values = solutions["apex"].vector().array()
    apex_values[d.dof_to_vertex_map(V)] = solutions["apex"].vector().array()
    ind_apex_max = apex_values.argmax()
    apex_coord = coords[ind_apex_max,:]

    # Update rhs
    L = v*d.Constant(0)*d.dx
    
    d.info("  Apex coord: ({0}, {1}, {2})".format(*apex_coord))
    d.info("  Num coords: {0}".format(len(coords)))
    d.info("  Num cells: {0}".format(len(cells_info)))
    
    # Calculate volume
    volume = 0.0
    for cell in d.cells(mesh):
        volume += cell.volume()
    
    d.info("  Volume: {0}".format(volume))
    d.info("")
    
    # Cases
    # =====
    #
    # 1) base: 1, apex: 0
    # 2) lv: 1, rv, epi: 0
    # 3) rv: 1, lv, epi: 0
    # 4) epi: 1, rv, lv: 0
    
    class ApexDomain(d.SubDomain):
        def inside(self, x, on_boundary):
            return d.near(x[0], apex_coord[0]) and d.near(x[1], apex_coord[1]) and \
                   d.near(x[2], apex_coord[2])
    
    apex_domain = ApexDomain()
    
    # Case 1:
    Poisson = 1
    DBC_11 = d.DirichletBC(V, 0, apex_domain, "pointwise")

    # Using Poisson
    if Poisson:
        d.solve(a==L, solutions["apex"], [DBC_10, DBC_11],
                solver_parameters={"linear_solver": "gmres"})

    # Using Eikonal equation
    else:
        Le = v*d.Constant(1)*d.dx
        d.solve(a==Le, solutions["apex"], DBC_11,
                solver_parameters={"linear_solver": "gmres"})

        # Create Eikonal problem
        eps = d.Constant(mesh.hmax()/25)
        y = solutions["apex"]
        F = d.sqrt(d.inner(d.grad(y), d.grad(y)))*v*d.dx - \
            d.Constant(1)*v*d.dx + eps*d.inner(d.grad(y), d.grad(v))*d.dx
        d.solve(F==0, y, DBC_11, solver_parameters={
            "linear_solver": "lu",
            "newton_solver": {"relative_tolerance":1e-5}})
    
    # Check that solution of the three last cases all sum to 1.
    sol = solutions["apex"].vector().copy()
    sol[:] = 0.0
    
    # Iterate over the three different cases
    for case in cases:

        # Solve linear system
        bcs = [d.DirichletBC(V, 1 if what == case else 0, \
                             facet_markers, markers[what], "topological") \
               for what in cases]
        d.solve(a==L, solutions[case], bcs, **solver_param)

        # Enforce bound on solution:
        solutions[case].vector()[solutions[case].vector()<d.DOLFIN_EPS] = d.DOLFIN_EPS
        solutions[case].vector()[solutions[case].vector()>1.0-d.DOLFIN_EPS] = 1.0-d.DOLFIN_EPS

        sol += solutions[case].vector()
    
    assert np.all(sol>.999), "Sum not 1..."

    # Return the solutions
    return solutions



def dolfin_to_vtk(output, basename):
    """
    Output a function to a basic vtk file
    """
    assert isinstance(output, d.Function)

    size = output.function_space().num_sub_spaces()
    if size>1:
        base_func = output.sub(0, True)
        base_space = base_func.function_space()
        assert base_space.num_sub_spaces() == 0, "expected only one subspace hiearchy"
        for i in range(1, size):
            assert base_space.ufl_element() == \
                   output.function_space().sub(i).ufl_element(), \
                   "expected sub spaces to be the same"

    else:
        base_space = output.function_space()
        
    lagrange = "Lagrange" in base_space.ufl_element().family()
        
    # If lagrange function space just save using dolfin pvd output
    if lagrange and base_space.ufl_element().degree()<=1:
        d.File("{}.pvd".format(basename)) << output
        return 

    coords = base_space.dofmap().tabulate_all_coordinates(base_space.mesh())
    coords = np.reshape(coords, (len(coords)/3, 3))
    dofs = output.vector().array()
    dofs = np.reshape(dofs, (len(dofs)/size, size))

    dof_format = " ".join("{}" for _ in range(size))+"\n"

    # write .vtk-file
    with open("{}.vtk".format(basename), "w") as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('VTK from dolfin\n')
        f.write('ASCII\n')
        f.write('DATASET UNSTRUCTURED_GRID\n')
        f.write('POINTS {} float\n'.format(len(coords)))
    
        for coord in coords:
            f.write('%f %f %f\n' % tuple(coord))
    
        f.write('POINT_DATA {}\n'.format(len(coords)) )
        f.write('VECTORS vectors float\n')
        for dof in dofs:
            f.write(dof_format.format(*tuple(dof)))
    
def dolfin_to_carpfile(mesh, basename, markers=None, \
                       vert_fields=None, cell_fields=None):
    """
    NOT DEBUGGED:
    Write carp mesh and fields to file from dolfin data

    mesh : dolfin.Mesh
        The dolfin.mesh which should be written to file
    basename : str
        Basename of file which all data will be written to
    markers : dict (optional)
        A dict of name to markers of facet booundaries contained in the mesh
    vert_fields : dict (optional)
        A dict between named vertex field data and dolfin Functions
    cell_fields : dict (optional)
        A dict between named cell field data and dolfin Functions
    
    """
    import dolfin as d
    import numpy as np

    boundary = d.CompiledSubDomain("on_boundary")

    d.warning("This function is not tested...")

    boundary_facets = d.FacetFunction("size_t", mesh, 0)
    boundary.mark(boundary_facets, 1)
    num_boundary_facets = np.sum(boundary_facets.array()==1)
    
    with open(basename+".pts", "w") as f:
        f.write("{0}\n".format(mesh.num_vertices()))
        [f.write("{0:.10f} {1:.10f} {2:.10f}\n".format(*coord)) \
         for coord in mesh.coordinates()]

    with open(basename+".elem", "w") as f:
        f.write("{0}\n".format(mesh.num_cells()))
        [f.write("Tt {0} {1} {2} {3} 0\n".format(*cell)) \
         for cell in mesh.cells()]
    
    with open(basename+".surf", "w") as f:
        f.write("{0}\n".format(num_boundary_facets))
        [f.write("Tr {0} {1} {2}\n".format(*facet.entities(0))) \
         for facet in d.SubsetIterator(boundary_facets, 1)]

    # If generating mapping between vertices and boundaries
    if markers:

        # Get the facet markers
        facet_markers = mesh.domains().facet_domains(mesh)

        # Iterate over markers
        for boundary_name, marker in markers.items():
            
            vertices = set()
            for face in SubsetIterator(facet_markers, marker):
                vertices.update(face.entities(0))
        
            with open(basename+"_"+boundary_name+".vtx", "w") as f:
                f.write("{0:d}\nintra \n".format(int(len(vertices))))
                [f.write("{0}\n".format(vert)) for vert in vertices]

        # Apex node...
        #with open(basename+"_apex.vtx", "w") as f:
        #    f.write("{0:d}\nintra \n".format(1))
        #    f.write("{0:d}\n".format((apex_point.array()==1).nonzero()[0][0]))

    # If outputing vertex fields
    if vert_fields:

        # Get dof mapping
        dofs_to_vert, vectordofs_to_vert, vectordofs_to_subvert = \
                      dofs_to_verts(mesh)

        # Iterate over the passed fields
        for field_name, func in vert_fields.items():
            values = func.vector().array()

            # If scalar field
            if mesh.num_vertices() == len(dofs):
                reordered_values = values[dofs_to_vert]

                # Write the field to file
                with open(basename+"_"+field_name+".dat", "w") as f:
                    [f.write("{0:.10f}\n".format(value)) for value in reordered_values]
            
            # If vector field
            elif mesh.num_vertices() == 3*len(dofs):
                raise NotImplementedError
            else:
                raise ValueError("Field and mesh do not match: "+field_name)
    
