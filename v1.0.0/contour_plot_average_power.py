from __future__ import division
import pylab as plt
from utils import *
import pandas as pd

df = pd.read_csv("max_min_coupling_v1.txt",delimiter=" ")
df.to_csv( "max_min_coupling_v1.csv", encoding='utf-8', index=False)          

WGLength, TLength, average_power = read_data(
    filename=r'max_min_coupling_v1',
    column=[0, 1, 3])

#SplittingEfficiency = Power_1/(Power_1 + Power_2)

num_pts = 101

WGLength = np.array(WGLength).reshape(num_pts, num_pts)
TLength = np.array(TLength).reshape(num_pts, num_pts)
average_power = np.array(average_power).reshape(num_pts, num_pts)*100

x, y = WGLength, TLength

fig0 = plt.figure()

levels = [0, 40, 50, 60, 100]
contour = plt.contour(x, y, average_power, levels, colors='k')
plt.clabel(contour, colors='k', fmt='%2.1f', fontsize=12)
contour_filled = plt.contourf(x, y, average_power, levels) #cmap='viridis')
plt.colorbar(contour_filled)

plt.xlabel(r'Waveguide Length [$\mu$m]')
plt.ylabel(r'Taper/Phase Compensation Section Length [$\mu$m]')
plt.title(r'The average power of port 1 in $\lambda$ 1500-1600 nm [%]')
fig0.savefig('figures/average_power.png', dpi=400)  # save the figure to file