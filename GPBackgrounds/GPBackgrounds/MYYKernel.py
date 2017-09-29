from sklearn.gaussian_process.kernels import *
from scipy.spatial.distance import squareform

def _convert_to_double(X):
    return np.ascontiguousarray(X, dtype=np.double)

class Gibbs(Kernel):
    """Gibbs kernel

    sqrt( (2*l(x)l(x')) / (l(x)**2 + l(x')**2 ) * exp( -(x-x')**2 / (l(x)**2 + l(x')**2) )
    where, l(x) = b*x + c

    This is the RBF kernel where the length scale is allowed to vary by some
    function, in this case a linear function.

    """

    def __init__(self, b=1.0, b_bounds=(1e-5,1e5), c=1.0, c_bounds=(1e-5,1e5)):
        self.b = b
        self.b_bounds = b_bounds
        self.c = c
        self.c_bounds = c_bounds

    @property
    def hyperparameter_b(self):
        return(Hyperparameter("b","numeric",self.b_bounds))

    @property
    def hyperparameter_c(self):
        return(Hyperparameter("c","numeric",self.c_bounds))



    def __call__(self, X, Y=None, eval_gradient=False):
        """Retrn K(X,Y)

        """
        def l(x):  #Helper function
            return self.b * x + self.c

        X = np.atleast_2d(X)

        s = X.shape
        if len(s) != 2:
            raise ValueError('A 2-dimensional array must be passed.')

        if Y is None:
            m, n = s
            K = np.zeros((m * (m - 1)) // 2, dtype=np.double)
            X = _convert_to_double(X)
            t = 0
            for i in xrange(0, m - 1):
                for j in xrange(i + 1, m):
                    xi_xj = X[i] - X[j]
                    li = l(X[i])  # b*x+c
                    lj = l(X[j])
                    li2_lj2 = np.dot(li,li) + np.dot(lj,lj)  # l(x)^2 + l(x')^2
                    coeff = np.sqrt(2*li*lj/(li2_lj2))
                    K[t] = coeff * np.exp(-1*(xi_xj*xi_xj) / li2_lj2 )
                    t = t + 1

            K = squareform(K)
            np.fill_diagonal(K,1)

            if eval_gradient:
                # approximate gradient numerically
                def f(theta):  # helper function
                    return self.clone_with_theta(theta)(X, Y)
                return K, _approx_fprime(self.theta, f, 1e-10)
            else:
                return K

        else:
            mx, nx = s
            sy = Y.shape
            my, ny = sy
            K = np.zeros((mx, my), dtype=np.double)
            X = _convert_to_double(X)
            Y = _convert_to_double(Y)
            t = 0
            for i in xrange(0, mx):
                for j in xrange(0, my):
                    xi_yj = X[i] - Y[j]
                    li = l(X[i])  # b*x+c
                    lj = l(Y[j])
                    li2_lj2 = li*li + lj*lj  # l(x)^2 + l(x')^2
                    coeff = np.sqrt(2*li*lj/(li2_lj2))
                    K[i][j] = coeff * np.exp(-1*(xi_yj* xi_yj) / li2_lj2 )
                    t = t + 1

            #K = squareform(K)
            #np.fill_diagonal(K,1)

            if eval_gradient:
                # approximate gradient numerically
                def f(theta):  # helper function
                    return self.clone_with_theta(theta)(X, Y)
                return K, _approx_fprime(self.theta, f, 1e-10)
            else:
                return K

        pass # __call__

    def diag(self, X):
        return np.diag(self(X))

    def is_stationary(self):
        return False

    def __repr__(self):
        return "{0}(b={1:.3g}, c={2:.3g})".format(
                self.__class__.__name__, self.b, self.c)


class FallExp(Kernel):
    """Falling exponential kernel

    exp( (d - (x+x'))/(2*a) )

    """
    def __init__(self, d=1.0, d_bounds=(1e-5,1e5), a=1.0, a_bounds=(1e-5,1e5)):
        self.d = d
        self.d_bounds = d_bounds
        self.a = a
        self.a_bounds = a_bounds

    @property
    def hyperparameter_d(self):
        return(Hyperparameter("d","numeric",self.d_bounds))

    @property
    def hyperparameter_a(self):
        return(Hyperparameter("a","numeric",self.a_bounds))


    def __call__(self, X, Y=None, eval_gradient=False):
        """Return K(X,Y)

        """
        X = np.atleast_2d(X)

        s = X.shape
        if len(s) != 2:
            raise ValueError('A 2-dimensional array must be passed.')

        if Y is None:
            m, n = s
            K = np.zeros((m * (m - 1)) // 2, dtype=np.double)
            X = _convert_to_double(X)
            t = 0
            for i in xrange(0, m - 1):
                for j in xrange(i + 1, m):
                    xi_xj = np.dtype('d')
                    xi_xj = X[i] + X[j]
                    K[t] = np.exp( (self.d - xi_xj) / (2*self.a) )
                    t = t + 1

            K = squareform(K)
            np.fill_diagonal(K,1)

            if eval_gradient:
                # approximate gradient numerically
                def f(theta):  # helper function
                    return self.clone_with_theta(theta)(X, Y)
                return K, _approx_fprime(self.theta, f, 1e-10)
                return K, None
            else:
                return K

        else:
            mx, nx = s
            sy = Y.shape
            my, ny = sy
            K = np.zeros((mx , my), dtype=np.double)
            X = _convert_to_double(X)
            Y = _convert_to_double(Y)
            for i in xrange(0, mx):
                for j in xrange(0, my):
                    xi_yj = np.dtype('d')
                    xi_yj = X[i] + Y[j]
                    K[i][j] = np.exp( (self.d - xi_yj) / (2*self.a) )

            #K = squareform(K)

            if eval_gradient:
                # approximate gradient numerically
                def f(theta):  # helper function
                    return self.clone_with_theta(theta)(X, Y)
                return K, _approx_fprime(self.theta, f, 1e-10)
            else:
                return K

        pass # __call__

    def diag(self,X):
        return np.diag(self(X))

    def is_stationary(self):
        return False

    def __repr__(self):
        return "{0}(d={1:.3g}, a={2:.3g})".format(
                self.__class__.__name__, self.d, self.a)



    pass # fallExp


# adapted from scipy/optimize/optimize.py for functions with 2d output
def _approx_fprime(xk, f, epsilon, args=()):
    f0 = f(*((xk,) + args))
    grad = np.zeros((f0.shape[0], f0.shape[1], len(xk)), float)
    ei = np.zeros((len(xk), ), float)
    for k in range(len(xk)):
        ei[k] = 1.0
        d = epsilon * ei
        grad[:, :, k] = (f(*((xk + d,) + args)) - f0) / d[k]
        ei[k] = 0.0
    return grad
