import numpy as np
import cheby as c
import scipy.integrate as integrate
import matplotlib.pyplot as plt

f = lambda x: np.cos(np.pi * x / 2)**3 + ((x + 1)**3)/8

def projection(N, M, n, x):
    """
    Parameters
    ----------
    N : Columns - total of N polynomials
    M : Rows - powers of x
    n : nth order of polynomial
    
    Returns
    -------
      <f|pn>
    ---------- x pn
      <pn|pn>
    """
    num = integrate.quad(lambda x: c.Cheby(N, M).eval_coeffs_pn(x, n) * f(x),
                                -1, 1)[0]
    den = integrate.quad(lambda x: c.Cheby(N, M).eval_coeffs_pn(x, n) * \
                                 c.Cheby(N, M).eval_coeffs_pn(x, n),
                                -1, 1)[0]
    fn_hat = num/den
    return fn_hat * c.Cheby(N, M).eval_coeffs_pn(x, n)

K = 4
x = np.linspace(-1, 1, 1000)
Pf = np.zeros(len(x))
for i in range(K):
    for j in x:
        Pf[i] += projection(8, 8, i, j)
print(f(x), Pf)
plt.plot(x, f(x), label="f(x)")
plt.plot(x, Pf, label="P$_N$f")
plt.legend()
plt.title("N=4")
plt.show()