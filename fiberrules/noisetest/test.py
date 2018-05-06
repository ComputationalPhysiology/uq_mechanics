import numpy as np
from dolfin import *


from kle import KLE
from covariance import CovarianceMatrix
from kernel import Matern





#mesh = UnitSquareMesh(10,10)
mesh =Mesh("LVentricle.xml")
#mesh.coordinates()[:] = 2.*mesh.coordinates()-1.

#Gaussian field with correlation length l
p , l = np.inf, 1
#number of terms in expansion
k = 20


kernel = Matern(p = p,l = l)	#0,1Exponential covariance kernel
kle  = KLE(mesh, kernel, verbose = True)
kle.compute_eigendecomposition(k = k)#20 by default

print "eigenvalues:", kle.l
plt.close('all')


#ploteigenvectors(kle)
#Generate realizations
noise = kle.realizations()


# The markers are stored within the mesh
#ffun = MeshFunction("size_t", mesh, 2, mesh.domains())
# The facet function contains some very large number,
# so if you want plot it you need to surpress those values
#ffun.array()[ffun.array() == max(ffun.array())] = 0
#plot(ffun,mesh, interactive =True))

#plotting the field
V = FunctionSpace(mesh, "CG", 1)
random_field =Function(V)
random_f = np.zeros((mesh.num_vertices()),dtype =float)
random_f = noise
random_field.vector()[:] = random_f[dof_to_vertex_map(V)]
plot(random_field, mesh, interactive =True)
