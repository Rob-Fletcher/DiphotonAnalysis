from sklearn.gaussian_process.kernels import *
from scipy.spacial.distance import squareform

def _convert_to_double(X):
    return np.ascontiguousarray(X, dtype=np.double)

class Gibbs(kernel):
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
        return(k.Hyperparamter("b","numeric",self.b_bounds))

    @property
    def hyperparameter_c(self):
        return(Hyperparamter("c","numeric",self.c_bounds))



    def __call__(self, X, Y=None, eval_gradient=False):
        """Retrn K(X,Y)

        """
        if not Y==None:
            raise("Y should be none.")

        def l(x):  #Helper function
            return self.b * x + self.c

        X = np.atleast_2d(X)

        s = X.shape
        if len(s) != 2:
            raise ValueError('A 2-dimensional array must be passed.')

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
                K[t] = coeff * np.exp(-1*np.dot(xi_xj, xi_xj) / li2_lj2 )
                t = t + 1

        K = squareform(K)
        np.fill_diagonal(K,1)

        if eval_gradient:
            raise("eval_gradient not implemented yet.")
        else:
            return K

        pass # __call__

    def __repr__(self):
        return "{0}(b={1:.3g}, c={2:.3g})".format(
                self.__class__.__name__, self.b, self.c)


class fallExp(kernel):
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
        return(k.Hyperparamter("d","numeric",self.d_bounds))

    @property
    def hyperparameter_a(self):
        return(Hyperparamter("a","numeric",self.a_bounds))



    def __call__(self, X, Y=None, eval_gradient=False):
        """Return K(X,Y)

        """
        X = np.atleast_2d(X)

        s = X.shape
        if len(s) != 2:
            raise ValueError('A 2-dimensional array must be passed.')

        m, n = s
        K = np.zeros((m * (m - 1)) // 2, dtype=np.double)
        X = _convert_to_double(X)
        t = 0
        for i in xrange(0, m - 1):
            for j in xrange(i + 1, m):
                xi_xj = X[i] + X[j]
                K[t] = np.exp( (self.d - xi_xj) / (2*self.a) )
                t = t + 1

        K = squareform(K)
        np.fill_diagonal(K,1)

        if eval_gradient:
            raise("eval_gradient not implemented yet.")
        else:
            return K

        pass # __call__

    def __repr__(self):
        return "{0}(d={1:.3g}, a={2:.3g})".format(
                self.__class__.__name__, self.d, self.a)



    pass # fallExp
