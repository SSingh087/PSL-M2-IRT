import numpy as np
import matplotlib.pyplot as plt 

xmin, xmax = -1, 1
tmin, tmax = 0, 5

x_all_cells, t_all_cells = 10, 10
X = np.linspace(xmin, xmax, x_all_cells, endpoint=True)
T = np.linspace(tmin, tmax, t_all_cells, endpoint=True)

RHO = np.zeros((len(T), len(X)))
V = np.zeros((len(T), len(X)))

delta_x = X[-1] - X[0] / len(X)

# declare initial array
for i in range(len(X)):
    if X[i] < 0:
        RHO[0, i] = 1 
        V[:, i] = 1
    elif X[i] > 0:
        RHO[0, i] = .5 
        V[:, i] = 0

def delta_t(i):
    return min(delta_x / max(V[i,],1e-5))

for i in range(0, len(T)-1):
    # boundary condition 
    RHO[i+1, 0] = RHO[i+1, 1]
    for j in range(1, len(X)-1):
        #calcualting v and rho at interface 
        #print(RHO[i, j])
        rho_n_l = RHO[i, j] - (RHO[i, j+1] - RHO[i, j-1]) / 4
        v_n_l = V[i, j] - (V[i, j+1] - V[i, j-1]) / 4
        rho_n_r = RHO[i, j] + (RHO[i, j+1] - RHO[i, j-1]) / 4
        v_n_r = V[i, j] + (V[i, j+1] - V[i, j-1]) / 4
        #updating values
        RHO[i+1, j] = RHO[i, j] - delta_t(i) / delta_x * (rho_n_r * v_n_r - rho_n_l * v_n_l)
    #print(delta_t(i))

for i in [0, 2, 4]:
    plt.plot(X, RHO[i, ], label=i)
    plt.legend()
plt.show()
