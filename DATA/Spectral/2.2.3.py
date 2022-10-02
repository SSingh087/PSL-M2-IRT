import numpy as np
import matplotlib.pyplot as plt
from waveletFunctions import wavelet
from scipy.io import readsav

data=readsav("mix_c1_20010331.sav")
t = data.get('t_sc')
t_hr = t/3600
bx, by, bz = data.get('bx'), data.get('by'), data.get('bz')

z = (t_hr >= 17.7) & (t_hr <=18.2)
t_hr = t[z]
Bx = bx[z]
dt = (t[-1] - t[0]) / len(t)

wave, period, scale, coi = wavelet(Bx, dt)
power = (np.abs(wave)) ** 2
print(min(scale), max(scale))
fig, ax = plt.subplots(2, figsize=(15,10))
ax[0].plot(t_hr/3600, Bx)
ax[0].set_ylabel('B$_x$ (nT)', fontsize=18)
ax[0].set_xlabel('time (hr)', fontsize=18)
ax[1].plot(t_hr/3600, coi, linewidth="2", color="black", linestyle="--")
ax[1].contourf(t_hr/3600, period, power**.1)
ax[1].set_yscale('log')
ax[1].set_ylabel('$\\tau$ (s)', fontsize=18)
ax[1].set_xlabel('time (hr)', fontsize=18)
ax[1].set_ylim(np.min(scale), np.max(coi))

plt.show()
