import numpy as np
import cheby as c
import scipy.integrate as integrate
import matplotlib.pyplot as plt

f = lambda x: np.cos(np.pi * x / 2)**3 + ((x + 1)**3)/8

def projection(N, n, x):
    """
    Parameters
    ----------
    N : Columns - total of N polynomials
    M : Rows - powers of x
    n : nth order of polynomial
    
    Returns
    -------
      <f|pn>
    ---------- * pn
      <pn|pn>
    """
    num = integrate.quad(lambda x: c.Cheby(N).eval_coeffs_pn(x, n) * f(x),
                                -1, 1)[0]
    den = integrate.quad(lambda x: c.Cheby(N).eval_coeffs_pn(x, n) * \
                                 c.Cheby(N).eval_coeffs_pn(x, n),
                                -1, 1)[0]
    fn_hat = num/den
    return fn_hat * c.Cheby(N).eval_coeffs_pn(x, n) / np.sqrt(1 - x**2)

N = 4
x = np.linspace(-1, 1, 100)
Pf = np.zeros(len(x))
for i in range(N+1):
    for j in x:
        Pf[i] += projection(N, i, j)
plt.plot(x, f(x), label="f(x)")
plt.plot(x, Pf, label="P$_N$f")
plt.xlabel('x')
plt.legend()
plt.title("N=4")
plt.show()
