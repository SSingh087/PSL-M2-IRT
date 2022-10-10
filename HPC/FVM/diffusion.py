import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
xmin, xmax = -1, 1
tmin, tmax = 0, 5

x_all_cells, t_all_cells = 100, 500
X = np.linspace(xmin, xmax, x_all_cells, endpoint=True)
T = np.linspace(tmin, tmax, t_all_cells, endpoint=True)
U = np.zeros((len(T), len(X)))
v = 0.1
#delta_t = T[-1] - T[0] / len(T)
delta_x = X[-1] - X[0] / len(X)
delta_t = lambda coef, v: coef * delta_x / v

val = v * (delta_t(0.8, 0.1) / delta_x)

# declare initial array
U[0, :] = np.exp(-X**2)
for i in range(len(T)-1):
    # boundary condition 
    U[i, 0] = U[i, 1]
    for j in range(1, len(X)):
        #updating values
        U[i+1, j] = U[i, j] - val * (U[i, j] - U[i, j-1])

for i in [0,20,50,70,100]:
    plt.plot(X, U[i, ], label=i)
    plt.legend()
#plt.show()
