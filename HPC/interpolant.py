import numpy as np
import cheby as c
import matplotlib.pyplot as plt

f = lambda x: np.cos(np.pi * x / 2)**3 + ((x + 1)**3)/8
           
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
    def __init__(self, N, M, K):
        self.N, self.M = N, M
        self.K = K
        self.wi = np.pi / K
        
    def gamma_n(self, n):
        """𝛾_𝑛 = 𝑁 ∑︁ 𝑝_𝑛(𝑥_𝑖)**2 * 𝑤_𝑖 
        for each collocation points for a given nth polynomial."""
        gn = 0
        for i in range(self.K+1):
            gn += c.Cheby(self.N, self.M).eval_coeffs_pn(xi(i, self.K), 
                                         n)**2 * self.wi
        return gn 
    
    def fn_tilde(self, n):
        """
        \tilde{𝑓_𝑛} = 1/𝛾_𝑛 ∑︁ 𝑓(𝑥_𝑖) * 𝑝_𝑛(𝑥_𝑖) * 𝑤_𝑖 
        for each collocation points for a given nth polynomial.
        """
        fn = 0
        for i in range(self.K+1):
            fn += f(xi(i, self.K)) * c.Cheby(self.N, 
                     self.M).eval_coeffs_pn(xi(i, self.K), n) * self.wi
        return fn/self.gamma_n(n)
        
    def do_interpolant(self, x):
        """
        ∑︁\tilde{𝑓_𝑛} * 𝑝_𝑛(𝑥)
        for each nth polynomial for a given x.
        """
        fn = np.zeros(self.K+1)
        pn = np.zeros(self.K+1)
        for i in range(self.K+1):
            fn[i] = self.fn_tilde(i) 
            pn[i] = c.Cheby(self.N, self.M).eval_coeffs_pn(x, i)
        return np.sum(fn*pn)

K = 8
N, M = 10, 10
x = np.linspace(-1,1,100)
inp = np.zeros(len(x))
for i in range(len(x)):
    inp[i] = interpolant(N, M, K).do_interpolant(x[i])
plt.plot(x,f(x), label="f(x)")
plt.plot(x,inp, label="P$_N$f")
plt.legend()
plt.title("N=8")
plt.show()

#eval_spectral convergence of I N f to f
