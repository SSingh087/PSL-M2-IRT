import numpy as np
import matplotlib.pyplot as plt
from waveletFunctions import wavelet
from scipy.interpolate import interp1d
from scipy.io import readsav

data=readsav("mix_c1_20010331.sav")
t = data.get('t_sc')
t_hr = t/3600
bx, by, bz = data.get('bx'), data.get('by'), data.get('bz')

z = (t_hr >= 17.7) & (t_hr <= 17.97)
t_hr = t[z]
Bx, By, Bz = bx[z], by[z], bz[z]

B = np.sqrt(bx**2 + by**2 + bz**2)[z]

FFT = lambda y: np.fft.fft(y)
Sn = lambda y: 2*T*abs(FFT(y))**2/N**2

T = (t_hr[-1] - t_hr[0])
N = len(t_hr)

f = (np.arange(int(N/2))+1)/T
print(max(f), min(f))
# 12.500797896107544 0.0016340912282493522

def interpolate(x, y):
    x_ = np.average(x)
    y_ = np.average(y)
    s1 = np.sum(x*y) - len(x)*x_*y_
    s2 = np.sum(x**2) - len(x)*x_**2

    a1 = s1/s2
    a0 = y_ - a1*x_
    print(a1)
    return a0+a1*x

#plt.plot(f, np.exp(interpolate(np.log(f), np.log(Sn(B)[:len(f)]))), color='blue')
#plt.plot(f, np.exp(interpolate(np.log(f), np.log(Sn(Bx)[:len(f)]))), linestyle="--", color='green')
#plt.plot(f, np.exp(interpolate(np.log(f), np.log(Sn(By)[:len(f)]))), linestyle="--", color='cyan')
#plt.plot(f, np.exp(interpolate(np.log(f), np.log(Sn(Bz)[:len(f)]))), linestyle="--", color='red')
plt.plot(f, Sn(B)[:len(f)], label="B$_{net}$", alpha=0.5, color='blue')
plt.plot(f, Sn(Bx)[:len(f)]+Sn(By)[:len(f)]+Sn(Bz)[:len(f)], label="B$_{sum}$", alpha=0.5, color='green')
plt.plot(f, np.exp(interpolate(np.log(f), np.log(Sn(B)[:len(f)]))), color='blue')
plt.plot(f, np.exp(interpolate(np.log(f), np.log(Sn(Bx)[:len(f)]+Sn(By)[:len(f)]+Sn(Bz)[:len(f)]))),
         linestyle="--", color='green')
#plt.plot(f, Sn(Bx)[:len(f)], linestyle="--", label="B$_{x}$", alpha=0.5, color='green')
#plt.plot(f, Sn(By)[:len(f)], linestyle="-.", label="B$_{y}$", alpha=0.5, color='cyan')
#plt.plot(f, Sn(Bz)[:len(f)], linestyle=":", label="B$_{z}$", alpha=0.5, color='red')
#plt.title("")"""
plt.xscale('log')
plt.yscale('log')
plt.xlabel('f (Hz)')
plt.ylabel('PSD')
plt.legend()
plt.show()
