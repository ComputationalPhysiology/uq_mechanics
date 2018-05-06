from fenics import *
from dolfin import *
from vtk import *
from tvtk.array_handler import *
from textwrap import dedent
import numpy as np
from scipy.spatial import KDTree
from distutils.version import StrictVersion

# Check for old vtk
vtk6 = StrictVersion(vtkVersion.GetVTKVersion()) >= StrictVersion("6.0.0")

def boundary_mesh(mesh, bfun, marker = None, u = None) :
    """
    Construct a new mesh from the boundary.
    """

    # boundary facets
    if marker :
        facets = np.where(bfun.array() == marker)[0]
    else :
        facets = np.where(bfun.array() != 0)[0]
    mesh.init(2, 0)
    conn = np.array(map(mesh.topology()(2, 0), facets))
    points = np.unique(conn.ravel())
    # renumber ids
    perm = { v: k for k, v in enumerate(points) }
    conn = np.vectorize(perm.get)(conn)
    coords = mesh.coordinates()[points, :]

    bmesh = Mesh(mesh.mpi_comm())
    dim = mesh.geometry().dim()
    editor = MeshEditor()
    editor.open(bmesh, dim-1, dim)

    editor.init_vertices(coords.shape[0])
    editor.init_cells(conn.shape[0])

    for vid, p in enumerate(coords) :
        editor.add_vertex(vid, p)
    for cid, c in enumerate(conn.view(np.uintp)) :
        editor.add_cell(cid, c)

    editor.close()

    return bmesh

def restrict_on_boundary(bmesh, u) :
    # new function on the boundary
    Vu = u.function_space()
    degree = Vu.ufl_element().degree()
    family = Vu.ufl_element().family()
    shape  = Vu.ufl_element().value_shape()
    if len(shape) == 0 :
        Vb = FunctionSpace(bmesh, family, degree)
    elif len(shape) == 1 :
        Vb = VectorFunctionSpace(bmesh, family, degree, shape[0])
    elif len(shape) == 2 :
        raise NotImplemented

    # brute-force
    gdim = bmesh.geometry().dim()
    dofb = Vb.dofmap().tabulate_all_coordinates(bmesh).reshape(-1, gdim)
    #dofb = Vb.tabulate_dof_coordinates().reshape(-1, gdim)#V.tabulate_dof_coordinates()
    dofu = Vu.dofmap().tabulate_all_coordinates(Vu.mesh()).reshape(-1, gdim)
    #dofu = Vu.tabulate_dof_coordinates().reshape(-1, gdim)

    if len(shape) > 0 :
        idxb = np.column_stack([ Vb.sub(i).dofmap().dofs() 
                    for i in xrange(0, shape[0]) ])
        idxu = np.column_stack([ Vu.sub(i).dofmap().dofs() 
                    for i in xrange(0, shape[0]) ])
    else :
        idxb = Vb.sub(i).dofmap().dofs()
        idxu = Vu.sub(i).dofmap().dofs()

    tree = KDTree(dofu[idxu[:, 0], :])
    src = tree.query(dofb[idxb[:, 0], :])[1]

    ub = Function(Vb)
    ub.vector()[np.array([1, 2])] = 1.0
    ub.vector()[:] = u.vector().array()[idxu][src, :].ravel()

    return ub

def dolfin2vtk(mesh, u = None) :

    # check for u in case of high order mesh
    gdim = mesh.geometry().dim()
    mdim = mesh.topology().dim()
    if isinstance(u, Function) :
        V = u.function_space()
        order = V.ufl_element().degree()
        # first we need to reorganize the entries of phi as matrix
        # whose column are components.
        # we do this by tabulating dofs of each subspace of V
        idx = np.column_stack([ V.sub(i).dofmap().dofs() 
                    for i in xrange(0, gdim) ])
        pts = V.dofmap().tabulate_all_coordinates(mesh).reshape(-1, gdim)
        #pts = V.tabulate_dof_coordinates().reshape(-1, gdim)
   	pts = pts[idx[:, 0], :] + u.vector().array()[idx]
        # always 3d coordinates
        pts = np.column_stack([pts, np.zeros((pts.shape[0], 3-gdim))]) 
        pts = array2vtkPoints(pts)

        # compute the connectivity from cell dofs
        reorder = np.array([[ 0, 1, 2 ],
                            [ 0, 1, 2, 5, 3, 4 ],
                            [ 0, 1, 2, 3, 9, 6, 8, 7, 5, 4 ]][mdim-1])
        if order == 1 :
            reorder = reorder[0:mdim+1]
        conn = np.array([ V.sub(0).dofmap().cell_dofs(c)[reorder]
                    for c in xrange(0, mesh.num_cells()) ])

        # now remap the ids with the new order
        perm = { v: k for k, v in enumerate(idx[:, 0]) }
        conn = np.vectorize(perm.get, otypes=[np.uintp])(conn)
        elms = array2vtkCellArray(conn)
    else :
        order = 1
        # coordinates
        pts = array2vtkPoints(mesh.coordinates())
        # connectivity
        dim = mesh.topology().dim()
        ndof = Cell(mesh, 0).num_vertices()
        conn = mesh.topology()(dim, 0)().reshape((-1, ndof))
        elms = array2vtkCellArray(conn)

    # grid
    grid = vtkUnstructuredGrid()
    grid.SetPoints(pts)

    # only these are supported by dolfin
    vtk_shape = { 1 : { 1 : VTK_LINE, 2 : VTK_TRIANGLE, 3 : VTK_TETRA },
                  2 : { 1 : VTK_QUADRATIC_EDGE,
                        2 : VTK_QUADRATIC_TRIANGLE,
                        3 : VTK_QUADRATIC_TETRA } }[order][mdim]
    grid.SetCells(vtk_shape, elms)

    return grid

def linear_grid_filter(grid, ndiv = 1) :
    # case of high-order
    geo = vtkUnstructuredGridGeometryFilter()
    if vtk6:
    	geo.SetInputData(grid)
    else:
    	geo.SetInput(grid)
        
    geo.Update()

    lingeo = vtkDataSetSurfaceFilter()
    lingeo.SetInputConnection(0, geo.GetOutputPort(0))
    lingeo.SetNonlinearSubdivisionLevel(1)
    lingeo.Update()

    return lingeo

def plot_vtk_grid(grid, ndiv = 0) :

    if ndiv > 0 :
        lingeo = linear_grid_filter(grid, ndiv)
        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(0, lingeo.GetOutputPort(0))
    else :
        mapper = vtkDataSetMapper()
        if vtk6:
           mapper.SetInputData(grid)
        else:
           mapper.SetInput(grid)

    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().EdgeVisibilityOn()

    renderer = vtkRenderer()
    renderer.AddActor(actor)

    renWin = vtkRenderWindow()
    renWin.AddRenderer(renderer)

    renWinInteractor = vtkRenderWindowInteractor()
    renWinInteractor.SetRenderWindow(renWin)

    style = vtk.vtkInteractorStyleTrackballCamera()
    renWinInteractor.SetInteractorStyle(style)

    renWin.Render()
    renWin.SetSize(800, 800)
    renWinInteractor.Start()

def surf_distance(geo1, geo2, surf = 'endo', u1 = None, u2 = None, plot = False) :
    
    #if not vtk6:
        #warning("Cannot compute surf distance. You need vtk version 5.9 or higher.")
    #    return -1.

    marker1 = { 'endo' : geo1.ENDO, 'epi' : geo1.EPI }[surf]
    marker2 = { 'endo' : geo2.ENDO, 'epi' : geo2.EPI }[surf]

    # restrict on boundary
    meshb1 = boundary_mesh(geo1.mesh, geo1.bfun, marker1)
    meshb2 = boundary_mesh(geo2.mesh, geo2.bfun, marker2)
    ub1 = restrict_on_boundary(meshb1, u1) if u1 else None
    ub2 = restrict_on_boundary(meshb2, u2) if u2 else None

    # vtk grid
    grid1 = dolfin2vtk(meshb1, ub1)
    grid2 = dolfin2vtk(meshb2, ub2)

    surface1_filter = vtkDataSetSurfaceFilter()
    surface1_filter.SetInputData(grid1)
    surface2_filter = vtkDataSetSurfaceFilter()
    surface2_filter.SetInputData(grid2)
    
    distance_filter = vtkDistancePolyDataFilter()
    distance_filter.SetInputConnection(0, surface1_filter.GetOutputPort())
    distance_filter.SetInputConnection(1, surface2_filter.GetOutputPort())
    distance_filter.ComputeSecondDistanceOn()
    distance_filter.SignedDistanceOff()
    distance_filter.Update()

    grid1_from_grid2 = distance_filter.GetOutput()\
            .GetPointData().GetScalars().GetRange()[1]
    grid2_from_grid1 = distance_filter.GetSecondDistanceOutput()\
            .GetPointData().GetScalars().GetRange()[1]

    dist = max(grid1_from_grid2, grid2_from_grid1)

    if plot :
        mapper1 = vtkPolyDataMapper()
        mapper1.SetInputConnection(distance_filter.GetOutputPort(0))
        distrange = distance_filter.GetOutput().GetPointData().GetScalars().GetRange()
        mapper1.SetScalarRange(distrange[0], distrange[1])
        mapper2 = vtkPolyDataMapper()
        mapper2.SetInputConnection(distance_filter.GetOutputPort(1))
        distrange = distance_filter.GetSecondDistanceOutput().GetPointData().GetScalars().GetRange()
        mapper2.SetScalarRange(distrange[0], distrange[1])

        actor1 = vtkActor()
        actor1.SetMapper(mapper1)
        actor1.GetProperty().SetRepresentationToWireframe()
        actor2 = vtkActor()
        actor2.SetMapper(mapper2)
        actor2.GetProperty().SetRepresentationToWireframe()

        renderer = vtkRenderer()
        renderer.AddActor(actor1)
        renderer.AddActor(actor2)

        renWin = vtkRenderWindow()
        renWin.AddRenderer(renderer)

        renWinInteractor = vtkRenderWindowInteractor()
        renWinInteractor.SetRenderWindow(renWin)

        renWin.Render()
        renWin.SetSize(400, 300)
        renWinInteractor.Start()

    return dist

def vtk_gradient_recovery(xyz, conn, u) :
    grid = vtk_create_grid(xyz, conn, u)

    gradient = vtkGradientFilter()
    if vtk6:
        gradient.SetInputData(grid)
    else:
        gradient.SetInput(grid)
    gradient.SetInputArrayToProcess(0, 0, 0, 0, 1)
    gradient.SetResultArrayName('grad')
    gradient.Update()
    grad_u = vtk2array(gradient.GetOutput().GetPointData().GetArray('grad'))

    return grad_u.reshape((-1, 3, 3))

def _plot_vtk_grid(grid) :
    mapper = vtkDataSetMapper()
    if vtk6:
      mapper.SetInputData(grid)
    else:
      mapper.SetInput(grid)

    actor = vtkActor()
    actor.SetMapper(mapper)
    #actor.GetProperty().SetRepresentationToWireframe()

    renderer = vtkRenderer()
    renderer.AddActor(actor)

    renWin = vtkRenderWindow()
    renWin.AddRenderer(renderer)

    renWinInteractor = vtkRenderWindowInteractor()
    renWinInteractor.SetRenderWindow(renWin)

    style = vtk.vtkInteractorStyleTrackballCamera()
    renWinInteractor.SetInteractorStyle(style)

    renWin.Render()
    renWin.SetSize(800, 800)
    renWinInteractor.Start()

def vtk_create_grid(xyz, conn, u) :
    # points
    pts = array2vtkPoints(xyz)
    
    # connectivity
    hexa_perm = np.array([ 0, 1, 3, 2, 4, 5, 7, 6 ])
    prism_perm = np.array([ 0, 3, 2, 4, 7, 6 ])
    isprism = conn[:,0] == conn[:,1]
    prism_conn = conn[isprism, :][:, prism_perm]
    hexa_conn = conn[np.logical_not(isprism), :][:, hexa_perm]
   
    elms = array2vtkCellArray([hexa_conn, prism_conn])
    vtk_type = [ VTK_HEXAHEDRON ] * hexa_conn.shape[0] \
             + [ VTK_WEDGE ] * prism_conn.shape[0]

    # grid
    grid = vtkUnstructuredGrid()
    grid.SetPoints(pts)
    grid.SetCells(vtk_type, elms)

    # scalar
    uvtk = array2vtk(u)
    if len(u.shape) > 1 :
        uvtk.SetName("dispacement")
        grid.GetPointData().SetVectors(uvtk)
    else :
        uvtk.SetName("stress")
        grid.GetPointData().SetScalars(uvtk)

    return grid

def export_vtk(fname, xyz, conn, u) :
    grid = vtk_create_grid(xyz, conn, u)
    writer = vtkXMLUnstructuredGridWriter()
    #if vtk6:
    writer.SetInputData(grid)
    #else:
    #writer.SetInput(grid)
    writer.SetFileName(fname)
    writer.Write()

def export_msh(fname, xyz, conn, **kwargs) :

    # gmsh counts from 1
    num_nodes = xyz.shape[0]

    nodes = '\n'.join('{} {} {} {}'.format(i+1, *p)
                      for i, p in enumerate(xyz))
    perm_hexa = np.array([ 0, 1, 3, 2, 4, 5, 7, 6 ])
    elements = '\n'.join('{} 5 2 99 2 {} {} {} {} {} {} {} {}'\
                         .format(i+1, *((c+1)[perm_hexa]))
                         for i, c in enumerate(conn) if c[0] != c[1])
    perm_prism = np.array([ 0, 3, 2, 4, 7, 6 ])
    elements += '\n'
    elements += '\n'.join('{} 6 2 99 2 {} {} {} {} {} {}'\
                          .format(i+1, *((c+1)[perm_prism]))
                          for i, c in enumerate(conn) if c[0] == c[1])
    
    content = dedent(\
    """\
    $MeshFormat
    2.2 0 8
    $EndMeshFormat
    $Nodes
    {num_nodes}
    {nodes}
    $EndNodes
    $Elements
    {num_elements}
    {elements}
    $EndElements\
    """).format(num_nodes = num_nodes, nodes = nodes,
                num_elements = conn.shape[0], elements = elements)

    for name, data in kwargs.iteritems() :
        num_comp = data.shape[1]
        values = '\n'.join(' '.join(['{}']*(num_comp+1)).format(i+1, *p)
                           for i, p in enumerate(data))
        content += dedent(\
        """
        $NodeData
        1
        "{name}"
        1
        0.0
        3
        0
        {num_comp}
        {num_nodes}
        {values}
        $EndNodeData
        """).format(name = name, num_nodes = num_nodes, 
                    num_comp = num_comp, values = values)

    with open(fname, "w") as f :
        f.write(content)

