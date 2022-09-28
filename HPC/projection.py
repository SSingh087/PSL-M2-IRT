import numpy as np
import cheby as c
import scipy.integrate as integrate

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
    ----------
      <pn|pn>
    """
    num = integrate.quad(lambda x: c.Cheby(N, M).eval_coeffs_pn(x, n) * f(x),
                                -1, 1)[0]
    den = integrate.quad(lambda x: c.Cheby(N, M).eval_coeffs_pn(x, n) * \
                                 c.Cheby(N, M).eval_coeffs_pn(x, n),
                                -1, 1)[0]
    fn_hat = num/den
    return fn_hat * c.Cheby(N, M).eval_coeffs_pn(x, n)

N = 4
x = np.linspace(0, 1, 100)
Pf = np.zeros(len(x))
for i in range(N):
    for j in x:
        Pf[i] += projection(5, 5, i, j)
print(Pf)