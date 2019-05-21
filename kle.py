from covariance import CovarianceMatrix
from dolfin import *
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import factorized, aslinearoperator, splu, eigsh
from time import time
import numpy as np
from petsc4py import PETSc
import scipy.sparse as sp_sparse

__all__ = ['KLE']


class KLE:
    """
	Class for generating eigenmodes of the Karhunen-Lo\`{e}ve 
		expansion corresponding to covariance kernel \kappa(\cdot,\cdot)

	Parameters:
	-----------
	mesh:	Mesh object generated from FEniCS
		Contains information about the mesh topology

	kernel: function
		A stationary, translation-invariant covariance kernel
		 $\kappa(x,y) = \kappa(||x-y||)$. See, covariance/kernel.py for examples

	verbose: bool, optional
		Displays timing information

	
	Attributes:
	-----------
	V:	FunctionSpace
		Linear Lagrange Finite Element Space. See dolfin

	Q:	CovarianceMatrix
		Implementation of the Covariance matrix class 
		{Dense Matrix, FFT and HMatrix}. Currently default is set to HMatrix
		FFT is restricted to regular grids		


	Methods:
	--------
	ComputeEigendecomposition()
	Realizations()


	Requires:
	---------
		
	dolfin, numpy, scipy, matplotlib


	References:
	-----------
	"""

    def __init__(self, mesh, kernel, verbose=False, params={}):
        self.mesh = mesh
        self.pts = mesh.coordinates()
        self.V = FunctionSpace(mesh, "Lagrange", 2)

        self.kernel = kernel
        self.verbose = verbose

        rkmax = 10 if 'rkmax' not in params else params['rkmax']
        eps = 1.e-12 if 'eps' not in params else params['eps']

        start = time()
        self.Q = CovarianceMatrix('Hmatrix', self.pts, kernel,
                                  rkmax=rkmax, eps=eps)
        if verbose:
            print("Time to setup H-matrix is %g" % (time() - start))

        start = time()
        self._BuildMassMatrix()
        if verbose:
            print("Time to setup Mass matrix is %g" % (time() - start))

    def _BuildMassMatrix2(self):
        """TODO: Docstring for getBmat. We are solving for CX = BX where C is the covariance matrix
        and B is just a mass matrix. Here we assemble B. This is a private function. DONT call this
        unless debuging.

        :returns: PETScMatrix B
        """

        # B matrix is just a mass matrix, can be easily assembled through fenics
        # however, the ordering in fenics is not the mesh ordering. so we build a temp matrix
        # then use the vertex to dof map to get the right ordering interms of our mesh nodes
        V = FunctionSpace(self.mesh, "CG", 1)
        # Define basis and bilinear form
        u = TrialFunction(V)
        v = TestFunction(V)
        a = u * v * dx
        B_temp = assemble(a)

        B = PETSc.Mat().create()
        B.setType('aij')
        B.setSizes(len(self.pts), len(self.pts))
        B.setUp()

        B_ij = B_temp.array()

        v_to_d_map = vertex_to_dof_map(V)

        print('---------------------------')
        print('---------------------------')
        print(' Building Mass Matrix ')
        print('---------------------------')
        print('---------------------------')
        for node_i in range(0, len(self.pts)):
            for node_j in range(node_i, len(self.pts)):
                B_ij_nodes = B_ij[v_to_d_map[node_i], v_to_d_map[node_j]]
                if B_ij_nodes > 0:
                    B.setValue(node_i, node_j, B_ij_nodes)
                    B.setValue(node_j, node_i, B_ij_nodes)

        B.assemblyBegin()
        B.assemblyEnd()
        print('---------------------------')
        print('---------------------------')
        print(' Finished Mass Matrix ')
        print('---------------------------')
        print('---------------------------')

        self.M = B
        return


    def _BuildMassMatrix(self):
        u = TrialFunction(self.V)
        v = TestFunction(self.V)

        m = u * v * dx
        ##M = uBLASSparseMatrix()
        M = PETScMatrix()
        assemble(m, tensor=M)

        # Convert into
        test = M.mat()
        row, col, data = test.getValuesCSR()

        ##(row, col, data) = M.data()   # get sparse data
        col = np.intc(col)
        row = np.intc(row)
        n = M.size(0)
        self.M = sp_sparse.csc_matrix((data, col, row), shape=(n, n), dtype='d')
        # M_scipy = sp_sparse.csc_matrix((data,col,row),shape=(n,n),dtype =np.float)
        # self.M = M_scipy
        return


    def _RandomizedGHEP(self, k, p=20, twopass=False):
        from eigen import RandomizedGHEP as GHEP
        verbose = self.verbose

        Qt = _MQM(self.Q, self.M)
        Mt = _Mass(self.M)
        qm = _QM(self.Q, self.M)

        l, v = GHEP(Qt, Mt, k=k, p=p, twopass=twopass, \
                    verbose=verbose, BinvA=qm, error=False)
        if verbose:
            print("Ax %d, Bx %d, B^{-1}x %d" %(Qt.mvcount, Mt.mvcount, Mt.scount))

        return l, v


    def _arpack(self, k):
        Qt = _MQM(self.Q, self.M)
        Mt = _Mass(self.M)
        Minv = _MassInv(self.M)

        verbose = self.verbose
        l, v = eigsh(aslinearoperator(Qt), M=aslinearoperator(Mt), \
                     k=k, which='LM', Minv=aslinearoperator(Minv))

        l = l[::-1]
        v = v[:, ::-1]

        if verbose:
            print("Ax %d, Bx %d, B^{-1}x %d" % (Qt.mvcount, Mt.mvcount, Minv.mvcount))
        return l, v


    def compute_eigendecomposition(self, k, method='Arpack', params={}):
        """
        Computes the k largest modes of the KLE

        Parameters:
        ------------
        k:	int
            Number of eigenmodes

        method:	string, optional
            Method to compute eigendecomposition.
            Options are either 'Arpack' [1,2] (default) or 'Randomized' [3]

        params:	dict, optional
            Parameters to control accuracy of randomized eigendecomposition
            keys: 'p' (int) is the oversampling factor and
            'twopass' (bool) determines whether to use
            two-pass (more accurate but expensive) or one-pass
            (less accurate and less expensive) [3].


        Returns:
        ---------

        l:	(k,) ndarray
            Contains the eigenvalues in descending order

        v:	(N,k) ndarray
            Contains the corresponding eigenvectors.
            N is the number of grid points.


        Raises:
        -------
        NotImplementedError
            When the right parameter for 'method' is not used

        References
            ----------
            .. [1] ARPACK Software, http://www.caam.rice.edu/software/ARPACK/
           .. [2] R. B. Lehoucq, D. C. Sorensen, and C. Yang,  ARPACK USERS GUIDE:
                   Solution of Large Scale Eigenvalue Problems by
             Implicitly Restarted Arnoldi Methods.
            SIAM, Philadelphia, PA, 1998.
        .. [3] A.K. Saibaba and P.K. Kitanidis:
            Randomized square-root free algorithms for
                generalized Hermitian eigenvalue problems.
            Preprint http://arxiv.org/abs/1307.6885

        """

        start = time()
        if method == 'Arpack':
            l, v = self._arpack(k)
        elif method == 'Randomized':

            twopass = 'True' if 'twopass' not in params else \
                params['twopass']
            p = 20 if 'p' not in params else params['p']
            l, v = self._RandomizedGHEP(k, p, twopass)

        else:

            raise NotImplementedError

        verbose = self.verbose
        if verbose:
            print("Time taken for computing eigendecomposition is %g" % (time() - start))

        self.l, self.v = l, v
        print(l)
        return(l, v)


    def realizations(self, eps):
        """
        Generates a realization from the KLE assuming a Gaussian process


        Returns:
        --------
        v:	(N,) ndarray
            Realization of the gaussian process

        """
        ##import chaospy as cp
        k = self.l.shape[0]
        ##joint = cp.J(cp.Normal(0,1),cp.Normal(0,1), cp.Normal(0,1),cp.Normal(0,1))
        # eps = joint.sample(1).T.reshape(4,)

        #  eps = np.random.randn(k)
        # eps = np.random.normal(0,9,k)

        # [x,y), rand*(y-x) + x
        # eps = np.random.randn(k)*(y-x) + x
        # eps = np.random.normal(x,y,k)

        l, v = self.l, self.v

        return np.sum(np.dot(v, np.diag(np.sqrt(l) * eps)), axis=1)


class _MQM:
    def __init__(self, Q, M):
        self.Q, self.M = Q, M
        n = M.shape[0]
        self.shape = (n, n)

        self.dtype = 'd'

        self.mvcount = 0
        self.scount = 0

    def matvec(self, x):
        M, Q = self.M, self.Q
        self.mvcount += 1
        return M * Q.matvec(M * x)

    def clear(self):
        self.mvcount = 0
        self.scount = 0
        return


class _QM:
    def __init__(self, Q, M):
        self.Q, self.M = Q, M
        n = M.shape[0]
        self.shape = (n, n)

        self.dtype = 'd'

        self.mvcount = 0
        self.scount = 0

    def matvec(self, x):
        M, Q = self.M, self.Q
        self.mvcount += 1
        return Q.matvec(M * x)

    def rmatvec(self, x):
        M, Q = self.M, self.Q
        self.mvcount += 1
        return M * Q.matvec(x)

    def clear(self):
        self.mvcount = 0
        self.scount = 0
        return


class _Mass:
    def __init__(self, M):
        self.M = M
        n = M.shape[0]
        self.shape = (n, n)
        self.fac = splu(M)
        self.dtype = 'd'

        self.mvcount = 0
        self.scount = 0

    def matvec(self, x):
        self.mvcount += 1
        return self.M * x

    def solve(self, x):
        self.scount += 1
        return self.fac.solve(x)

    def clear(self):
        self.mvcount = 0
        self.scount = 0
        return


class _MassInv(_Mass):
    def __init__(self, M):
        _Mass.__init__(self, M)

    def matvec(self, x):
        self.mvcount += 1
        return self.fac.solve(x)

    def solve(self, x):
        self.scount += 1
        return self.M * x


if __name__ == '__main__':
    n = 50
