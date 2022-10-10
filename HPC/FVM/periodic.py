import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
xmin, xmax = 0, 2*np.pi
tmin, tmax = 0, 5

c = 1
x_all_cells, t_all_cells = 100, 100
delta_x = (xmax-xmin) /(x_all_cells-2)
X = np.linspace(xmin-delta_x, xmax+delta_x, x_all_cells, endpoint=True)
T = np.linspace(tmin, tmax, t_all_cells, endpoint=True)
U = np.zeros((len(T), len(X)))

delta_t = T[-1] - T[0] / len(T)
delta_x = X[-1] - X[0] / len(X)
val = c * (delta_t / delta_x)

# declare initial array
U[0, :] = np.sin(X)
for i in range(len(T)-1):
    # boundary condition 
    U[i, 0] = U[i, -2]
    for j in range(1, len(X)):
        #updating values
        U[i+1, j] = U[i, j] - val * (U[i, j] - U[i, j-1])

plt.plot(X, U[4, ])