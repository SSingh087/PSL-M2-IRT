import numpy as np
#import matplotlib.pyplot as plt
import scipy.integrate as integrate


K=4 # N in literature
#w0, wn = np.pi/2*K, np.pi/2*K
wi = np.pi/K
class Cheby():
    """
    Some used Parameters
    ----------
    x : point on x axis at which value has to be 
        calculated
    n : nth polynomial
    """
    def __init__(self, N, M):
        """
        Parameters
        ----------
        N : Columns - total of N polynomials
        M : Rows - powers of x
        """
        self.T = np.zeros((N, M))
        self.T[0,0] = 1
        self.T[1,1] = 1
    
    def eval_matrix(self):
        """
        returns the matrix $T_n$
        """
        for i in range(2,N):
            for j in range(M):
                self.T[i,j] = 2*self.T[i-1,j-1] - self.T[i-2,j]
        return self.T

    def eval_coeff(self, x, n):  
        """
        Returns
        -------
        sum of coefficients
        """
        coeff = 0
        if n > M:
            raise Exception("coefficients are not defined for \
                            nth power")
        for j in range(M):
            coeff += self.T[n,j]*x**j
        return(coeff)
        
    def get_func(self, x, func, n):
        """
        Parameters
        ----------
        func : function with which inner product is taken

        Returns
        -------
        inner product <pn|fn>
        """
        return func * self.eval_coeff(x, n)
    
    def self_prod(self, x, n):
        """
        Parameters
        ----------
        x : point on x axis at which value has to be 
            calculated
        n : nth polynomial

        Returns
        -------
        inner product with self <pn|pn>
        """
        return self.eval_coeff(x, n) * self.eval_coeff(x, n)
    
    def projection_coeff(self, n):
        """
        Returns
        -------
        f_hat = <pn|fn> / <pn|pn>
        """
        num = integrate.quad(self.get_func, -1, 1, args=(f, n))[0]
        den = integrate.quad(self.self_prod, -1, 1, args=(n))[0]
        f_hat = num/den
        return f_hat
    
    def projection(self, x, n):
        """
        Returns
        -------
        f_hat = <pn|fn> / <pn|pn> for a given value of x

        """
        return self.projection_coeff(n) * self.eval_coeff(x, n)
        
    def integrant(self):
        None
        
       
f = lambda x: np.cos(np.pi * x / 2)**3 + ((x + 1)**3)/8
           
def orth_check(x):
    """
    to check the orthogonality of the polynomials
    wrt to measure
    Parameters
    ----------
    x : point on x axis at which value has to be 
        calculated 
    Returns
    -------
    inner product.

    """
    n, m = 1, 4  #order of polynomial
    N = 8
    M = 10
    c = Cheby(N, M)
    Tn = c.eval_coeff(x, n)
    Tm = c.eval_coeff(x, m)
    return Tn * Tm / np.sqrt(1 - x**2)


if __name__ == '__main__':
    
    N = 8
    M = 10
    K = 4 
    cheby = Cheby(N, M)
    #print(cheby.eval_matrix())
    res, error = integrate.quad(orth_check, -1, 1)
    #print(res)
    Pn = 0
    x = 2
    for i in range(K):
        Pn += cheby.projection(x, i)
    print(Pn)

a = """gamma = 0
for i in range(K):
	xi = np.cos(np.pi * i / K)
	gamma += np.sum(compute(xi, 0)**2 * wi)
#for i in range(
#fn = 1/gamma * np.sum()
#print(gamma)

#print(T[7])"""


