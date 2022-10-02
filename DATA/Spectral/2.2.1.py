from scipy.io import readsav
import matplotlib.pyplot as plt
import numpy as np

data = readsav("mix_c1_20010331.sav")
t = data.get('t_sc')
bx, by, bz = data.get('bx'), data.get('by'), data.get('bz')
#print(t[155000]/3600)
#print(t[199000]/3600)
#print(t[180000]/3600)
plt.figsize=(15,10)
plt.plot(t/3600, bx, alpha=.8, label='B$_x$')
plt.plot(t/3600, by, alpha=.8, label='B$_y$')
plt.plot(t/3600, bz, alpha=.8, label='B$_z$')
plt.tight_layout()
plt.legend(ncol=3)
plt.xlabel('t (hr)', fontsize=18)
plt.ylabel('B (nT)', fontsize=18)
plt.xlim(17.7, 18.2)
plt.show()

B = np.sqrt(bx**2 + by**2 + bz**2)
plt.plot(t/3600, B)
plt.xlabel('t (hr)', fontsize=18)
plt.ylabel('B$_{net}$ (nT)', fontsize=18)
plt.xlim(17.7, 18.2)
plt.show()