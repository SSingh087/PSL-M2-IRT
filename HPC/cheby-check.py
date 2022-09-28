import cheby as c
import scipy.integrate as integrate

def cheby_matrix(N): 
    print(c.Cheby(N).eval_matrix())
    
def test_orthogonality(N, n, m):
    res, err = integrate.quad(c.Cheby(N).orth_func, -1, 1, args=(n, m))
    print(res)
    
def eval_coeffs(N, n, x):
    print(c.Cheby(N).eval_coeffs_pn(x, n))
    print(c.Cheby(N).eval_coeffs_p(x))
    
cheby_matrix(8)
eval_coeffs(8, 5, 2)
test_orthogonality(8, 1, 2)
test_orthogonality(8, 2, 2)