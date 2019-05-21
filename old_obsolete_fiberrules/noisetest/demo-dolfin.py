from dolfin import *
from fiberrules import *


parameters.allow_extrapolation = True

mesh = Mesh("LVentricle.xml")

fiber_space = FunctionSpace(mesh, "CG", 1)
fibers, sheet_normals, cross_sheet = \
        dolfin_fiberrules(mesh, fiber_space)

