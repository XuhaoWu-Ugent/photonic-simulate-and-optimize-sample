from __future__ import division
import pylab as plt
from utils import *
import pandas as pd


WGLength, TLength, Power_1, Power_2 = read_data(
    filename=r'BDC_pyfimm_v6',
    column=[0, 1, 3, 4])
#this "read_data" function is unnecessary. if for further commericial use,
#it could be replaced by: dataframe(data in pandas) to array: df = df.value
#array to dataframe: df = pd.DataFrame(df) while pandas is needed.

SplittingEfficiency = Power_1/(Power_1 + Power_2)

num_pts = 101 # need to be modified for different resolution of the data

WGLength = np.array(WGLength).reshape(num_pts, num_pts)
TLength = np.array(TLength).reshape(num_pts, num_pts)
SplittingEfficiency = np.array(SplittingEfficiency).reshape(num_pts, num_pts)*100

x, y = WGLength, TLength

fig0 = plt.figure()

levels = [0, 40, 50, 60, 100]
contour = plt.contour(x, y, SplittingEfficiency, levels, colors='k')
plt.clabel(contour, colors='k', fmt='%2.1f', fontsize=12)
contour_filled = plt.contourf(x, y, SplittingEfficiency, levels)#,cmap='viridis')
plt.colorbar(contour_filled)

plt.xlabel(r'Waveguide Length [$\mu$m]')
plt.ylabel(r'Taper/Phase Compensation Section Length [$\mu$m]')
plt.title(r'Splitting Efficiency at $\lambda$ = 1600 nm [%]')
fig0.savefig('figures/SplittingEfficiency @ 1600nm.png', dpi=400)  # save the figure to file

