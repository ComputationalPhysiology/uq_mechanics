from dolfin import *
from fiberrules import *

parameters.allow_extrapolation = True

mesh = Mesh("heart.xml.gz")

fiber_space = FunctionSpace(mesh, "Quadrature", 4)
fiber_space = FunctionSpace(mesh, "CG", 2)
fibers, sheet_normals, cross_sheet = \
        dolfin_fiberrules(mesh, fiber_space)

dolfin_to_vtk(fibers, "heart")
