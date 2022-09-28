import cheby as c
import scipy.integrate as integrate

def cheby_matrix(N, M): 
    print(c.Cheby(N, M).eval_matrix())
    
def test_orthogonality(N, M, x, n, m):
    res, err = integrate.quad(c.Cheby(N, M).orth_func, -1, 1, args=(n, m))
    print(res, err)
    
def eval_coeffs(N, M, n, x):
    print(c.Cheby(N, M).eval_coeffs_pn(x, n))
    print(c.Cheby(N, M).eval_coeffs_p(x))
    
#eval_coeffs(8, 8, 5, 1)
cheby_matrix(8, 10)
eval_coeffs(8, 10, 5, 2)