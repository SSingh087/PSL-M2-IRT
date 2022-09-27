import numpy as np
import matplotlib.pyplot as plt
from waveletFunctions import wavelet

delta_t, N = 1, 1000
T = N*delta_t

print("Min freq",1/T)
print("Max freq",N/(2*T))

t, n = np.arange(0, N, delta_t), np.arange(0, N/2)
f0 = 50/T
f1 = 30/T
print("Freq of signal",f0, f1)
y = 5*np.sin(2*np.pi*f0*t) + 4*np.sin(2*np.pi*f1*t)

dt = delta_t/T
print("dt:",dt)

wave, period, scale, coi = wavelet(y, dt)
power = (np.abs(wave)) ** 2


def DFT(n):
    un = 0
    for i in range(0,N-1):
        un += (1.0/N)*y[i]*np.exp(-2j*np.pi*i*n/N)
    return un

def PSD(n):
    return 2*T*abs(DFT(n))**2


fig, ax = plt.subplots(3, figsize=(15,8))
ax[0].plot(t, y, marker='+')
ax[0].set_ylabel('Amplitude', fontsize=18)
ax[0].set_xlabel('time (s)', fontsize=18)
ax[1].plot(n/T, PSD(n), marker='*')
ax[1].set_ylabel('PSD', fontsize=18)
ax[1].set_xlabel('freq (Hz)', fontsize=18)
ax[1].set_xlim(.02,.08)
ax[2].contourf(t, period, power)
ax[2].set_yscale('log')
ax[2].set_ylabel('freq (Hz)', fontsize=18)
ax[2].set_xlabel('time (s)', fontsize=18)
plt.show()
