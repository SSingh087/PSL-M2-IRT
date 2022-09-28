import numpy as np

class Cheby():
    """
    Some used Parameters
    ----------
    x : point on x axis at which value has to be 
        calculated
    n : nth order of polynomial
    """
    def __init__(self, N):
        """
        Parameters
        ----------
        N : Columns - total of N polynomials
        M : Rows - powers of x
        """
        self.N, self.M = N+1, N+1
        self.T = np.zeros((self.N, self.M))
        self.T[0,0] = 1
        self.T[1,1] = 1
        for i in range(2, self.N):
            for j in range(self.M):
                self.T[i,j] = 2*self.T[i-1,j-1] - self.T[i-2,j]

    def eval_matrix(self):
        """
        returns the matrix $T_n$
        """
        return self.T

    def eval_coeffs_pn(self, x, n):  
        """
        Returns
        -------
        sum of coefficients for given n
        """
        coeff = 0
        for j in range(self.M):
            coeff += self.T[n,j]*x**j
        return(coeff)
    
    def eval_coeffs_p(self, x):  
        """
        Returns
        -------
        sum of coefficients for all n
        """
        coeff = np.zeros(self.N)
        for i in range(self.N):
            for j in range(self.M):
                coeff[i] += self.T[i,j]*x**j
        return(coeff)
    
    def orth_func(self, x, n, m):
        """
        to check the orthogonality of the polynomials
        wrt to measure
        Parameters
        ----------
        x : point on x axis at which value has to be 
            calculated 
        n, m : order of polynomial
        Returns
        -------
        inner product.
    
        """
        Tn = self.eval_coeffs_pn(x, n)
        Tm = self.eval_coeffs_pn(x, m)
        return Tn * Tm / np.sqrt(1 - x**2)