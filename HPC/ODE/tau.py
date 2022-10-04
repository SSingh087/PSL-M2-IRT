import numpy as np
import interpolant as interpolant
import matplotlib.pyplot as plt
import cheby as c

source = lambda x: np.exp(x) - 4 * np.exp(1) / (1 + np.exp(1)**2)
analytic = lambda x: np.exp(x) - np.sinh(1.) / np.sinh(2.) * np.exp(2 * x) - \
                np.exp(1) / (1 + np.exp(1)**2)
der_func = lambda a, N : der(a, N, 2) - 4 * der(a, N, 1) + 4 * a

def der(a, N, df_order):
    b = np.zeros(N)
    if df_order == 1:
        for n in range(N):
            for p in range(n+1, N, 2):
                b[n] += p * a[p]
            if n!=0:
                b[n]*=2
    elif df_order == 2:
        for n in range(N):
            for p in range(n+2, N, 2):
                b[n] += p * (p**2 - n**2) * a[p]
            if n==0:
                b[n]*=.5
    return b

N = int(input("enter order: "))

Lij = np.zeros((N, N))
for i in range(N):
	a = np.zeros(N)
	a[i] = 1.
	Lij[:,i] = der_func(a, N)
print(Lij)

fns, fna = np.zeros(N), np.zeros(N)
for i in range(N):
    fns[i] = interpolant.interpolant(N).fn_tilde(i, source)
    fna[i] = interpolant.interpolant(N).fn_tilde(i, analytic)
    
for j in range(N):
	if (j%2 == 0): Lij[N-1, j] = 1
	else : Lij[N-2, j] = -1
fns[N-2] = 0

for j in range(N):
    Lij[N-1, j] = 1
fns[N-1] = 0

Lij_inv = np.linalg.solve(Lij, fns)

x = np.linspace(-1,1,100)
y = np.zeros(len(x))

T = c.Cheby(N)
for i in range(len(x)):
    res = 0
    for j in range(N):
        res += Lij_inv[j] * c.Cheby(N).eval_coeffs_pn(x[i], j)
    y[i] = res
plt.plot(x, y, label="tau method")
plt.plot(x, analytic(x), label="analytic")
plt.legend()
plt.xlim(-1,1)
plt.ylim(0, 0.5)
plt.title("N=12")
plt.show()

err = max(abs(y - analytic(x)))
print (N, " ", err) 