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
    
    def fn_tilde(self, n):
        """
        \tilde{𝑓_𝑛} = 1/𝛾_𝑛 ∑︁ 𝑓(𝑥_𝑖) * 𝑝_𝑛(𝑥_𝑖) * 𝑤_𝑖 
        for each collocation points for a given nth polynomial.
        """
        fn = 0
        for i in range(self.N+1):
            fn += f(xi(i, self.N)) * c.Cheby(self.N).eval_coeffs_pn(xi(i, 
                   self.N), n) * self.wi[i]
        return fn/self.gamma_n(n)
        
    def do_interpolant(self, x):
        """
        ∑︁\tilde{𝑓_𝑛} * 𝑝_𝑛(𝑥)
        for each nth polynomial for a given x.
        """
        fn = np.zeros(self.N+1)
        pn = np.zeros(self.N+1)
        for i in range(self.N+1):
            fn[i] = self.fn_tilde(i) 
            pn[i] = c.Cheby(self.N).eval_coeffs_pn(x, i)
        return np.sum(fn*pn)

N = 8
x = np.linspace(-1,1,100)
inp = np.zeros(len(x))
for i in range(len(x)):
    inp[i] = interpolant(N).do_interpolant(x[i])
plt.plot(x,f(x), label="f(x)")
plt.plot(x,inp, label="P$_N$f")
plt.legend()
plt.xlabel('x')
plt.title("N=8")
plt.show()
spectral_conv = abs(inp-f(x))
plt.title("N=8")
plt.ylabel("I$_N$f - f")
plt.xlabel('x')
plt.plot(x, spectral_conv)
plt.show()

