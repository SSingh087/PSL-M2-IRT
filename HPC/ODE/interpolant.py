import numpy as np
import cheby as c
           
xi = lambda i, K: np.cos(np.pi * i / K)

class interpolant():
    """
    Some used Parameters
    ----------
    x : point on x axis at which value has to be 
        calculated
    n : nth order of polynomial
    N : Columns - total of N polynomials
    M : Rows - powers of x
    """
    def __init__(self, N):
        self.N = N
        self.wi = np.zeros(self.N+1)
        for i in range(self.N+1):
            self.wi[i] = np.pi/self.N
        self.wi[0] /= 2
        self.wi[self.N] /= 2
        
    def gamma_n(self, n):
        """𝛾_𝑛 = 𝑁 ∑︁ 𝑝_𝑛(𝑥_𝑖)**2 * 𝑤_𝑖 
        for each collocation points for a given nth polynomial."""
        gn = 0
        for i in range(self.N+1):
            gn += c.Cheby(self.N).eval_coeffs_pn(xi(i, self.N), 
                                         n)**2 * self.wi[i]
        return gn 
    
    def fn_tilde(self, n, f):
        """
        \tilde{𝑓_𝑛} = 1/𝛾_𝑛 ∑︁ 𝑓(𝑥_𝑖) * 𝑝_𝑛(𝑥_𝑖) * 𝑤_𝑖 
        for each collocation points for a given nth polynomial.
        """
        fn = 0
        for i in range(self.N+1):
            fn += f(xi(i, self.N)) * c.Cheby(self.N).eval_coeffs_pn(xi(i, 
                   self.N), n) * self.wi[i]
        return fn/self.gamma_n(n)
    
    def do_interpolant(self, x, f):
        """
        ∑︁\tilde{𝑓_𝑛} * 𝑝_𝑛(𝑥)
        for each nth polynomial for a given x.
        """
        fn = np.zeros(self.N+1)
        pn = np.zeros(self.N+1)
        for i in range(self.N+1):
            fn[i] = self.fn_tilde(i, f) 
            pn[i] = c.Cheby(self.N).eval_coeffs_pn(x, i)
        return np.sum(fn*pn)