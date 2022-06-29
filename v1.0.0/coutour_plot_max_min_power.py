from __future__ import division
import pylab as plt
from utils import *
import pandas as pd

df = pd.read_csv("max_min_coupling.txt",delimiter=" ")
df.to_csv( "max_min_coupling.csv", encoding='utf-8', index=False)          

WGLength, TLength, min_max_coupling = read_data(
    filename=r'max_min_coupling',
    column=[0, 1, 2])

#SplittingEfficiency = Power_1/(Power_1 + Power_2)

num_pts = 101

WGLength = np.array(WGLength).reshape(num_pts, num_pts)
TLength = np.array(TLength).reshape(num_pts, num_pts)
max_min_coupling = np.array(min_max_coupling).reshape(num_pts, num_pts)

x, y = WGLength, TLength

fig0 = plt.figure()

levels = [0.1, 1, 20, 40]
contour = plt.contour(x, y, max_min_coupling, levels, colors='k')
plt.clabel(contour, colors='k', fmt='%2.1f', fontsize=12)
contour_filled = plt.contourf(x, y, max_min_coupling, levels) #cmap='viridis')
plt.colorbar(contour_filled)

plt.xlabel(r'Waveguide Length [$\mu$m]')
plt.ylabel(r'Taper/Phase Compensation Section Length [$\mu$m]')
plt.title(r'The ratio between max and min coupling in $\lambda$ 1500-1600 nm [dB]')
fig0.savefig('figures/min_max_coupling.png', dpi=400)  # save the figure to file