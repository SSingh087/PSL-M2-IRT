import numpy as np
import cheby as c
import matplotlib.pyplot as plt

f = lambda x: np.cos(np.pi * x / 2)**3 + ((x + 1)**3)/8
           
xi = lambda i, K: np.cos(np.pi * i / K)

class interpolant():
    def __init__(self, N, M, K):
        self.N, self.M = N, M
        self.K = K
        self.wi = np.pi / K
        
    def gamma_n(self, n):
        gn = 0
        for i in range(self.K):
            gn += c.Cheby(self.N, self.M).eval_coeffs_pn(xi(i, self.K), 
                                         n)**2 * self.wi
        return gn 
    
    def fn_tilde(self, n):
        fn = 0
        for i in range(self.K):
            fn += f(xi(i, self.K)) * c.Cheby(self.N, 
                     self.M).eval_coeffs_pn(xi(i, self.K), n) * self.wi
        return fn/self.gamma_n(n)
        
    def interpolant(self, x):
        I = 0
        for i in range(self.K):
            I += self.fn_tilde(i) * c.Cheby(self.N, 
                     self.M).eval_coeffs_pn(x, i)
        return I 



K = 8
N, M = 10, 10
x = np.linspace(-1,1,100)
inp = np.zeros(len(x))
for i in range(len(x)):
    inp[i] = interpolant(N, M, K).interpolant(x[i])
print(inp)
plt.plot(x,f(x), label="f(x)")
plt.plot(x,inp, label="P$_N$f")
plt.legend()
plt.title("N=4")
plt.show()
#gn = np.zeros(K)
#for i in range(K):
#    gn[i] = gamma_n(8, 10, 4, i)
#print(gn)
    