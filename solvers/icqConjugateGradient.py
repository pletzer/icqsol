from __future__ import print_function
import numpy


class ConjugateGradient:

    def __init__(self, mat, b):
        """
        Constructor
        @param mat dense, square matrix
        @param b right hand side vector
        """
        self.mat = mat
        self.b = b
        n = len(b)
        self.maxNumIters = n
        self.tol = 1.e-10
        self.verbose = False
        self.precond = numpy.array([mat[i, i] for i in range(n)])

    def setTolerance(self, tol):
        """
        Set tolerance
        @param tol tolerance
        """
        self.tol = tol

    def setMaxNumberOfIterations(self, maxNumIters):
        """
        Set maximum number of iterations
        @param maxNumIters  number of iterations
        """
        self.maxNumIters = maxNumIters

    def setVerbosity(self, verbose):
        """
        Set verbosity
        @param verbose True will print(out messages)
        """
        self.verbose = verbose

    def setDiagonalPreconditioner(self, precond):
        """
        Set the diagonal preconditioner
        @param precond preconditioner
        """
        self.precond = precond

    def solve(self, x0):
        """
        Solve linear system
        @param x0 initial guess for solution
        @return solution, error, and number of iterations
        """

        x = x0.copy()
        r = self.b - self.mat.dot(x)
        w = r / self.precond
        p = numpy.zeros(x0.shape, numpy.float64)
        beta = 0.0
        rho = r.dot(w)
        err = numpy.linalg.norm(r)
        k = 0
        while abs(err) > self.tol and k < self.maxNumIters:
            p = w + beta*p
            z = self.mat.dot(p)
            alpha = rho / p.dot(z)
            r -= alpha*z
            w = r / self.precond
            rhoOld = rho
            rho = r.dot(w)
            x += alpha*p
            beta = rho/rhoOld
            err = numpy.linalg.norm(self.b - self.mat.dot(x))
            if self.verbose:
                print('iteration {0} error = {1}'.format(k, err))
            k += 1

        return x, err, k

    def getSolutionError(self, vec):
        """
        Get the solution error
        @param vec solution
        @return error
        """
        return numpy.linalg.norm(self.b - self.mat.dot(vec))
