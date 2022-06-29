import subprocess
from time import gmtime, strftime, sleep

from scipy.interpolate import CubicSpline

import pyfimm as pf  # Every script must begin with this line
from utils import *
import numpy as np
subprocess.Popen(
    r'"C:\Program Files (x86)\PhotonD6.6.2\Fimmwave\bin64\fimmwave.exe" -pt 5101')  # Open fimmwave software.
sleep(10)  # Wait for a couple of minutes to let App open.
print "Start"
pf.connect()  # this connects to the FimmWave application.
begin_time = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
print "This program starts at {}".format(begin_time)

# Set Parameters (Your copy of FIMMWAVE has default values for these. You can change more than shown here. See __jaredwave.py
import os

ScriptPath, ScriptFile = os.path.split(os.path.realpath(__file__))  # Get directory of this script

# widths
WG_width = 0.3
WG_gap = 0.2
Clad_width = 1

# Read data from the optimization results
df = pd.read_csv("test_results_xuhao.txt", delimiter=" ")
df.to_csv("test_results_xuhao.csv", encoding='utf-8', index=False)
Opt_results = pd.read_csv("test_results_xuhao.csv")
# Opt_results_ideal = Opt_results[(Opt_results["result"] >0.49) & (Opt_results["result"] < 0.51)] #for 100/0 BDC
Opt_results_ideal = Opt_results[Opt_results["result"] < 0.04]
Opt_results_ideal.to_csv("Opt_results_ideal.csv", encoding='utf-8', index=False)
wg_length, T_length, y_upper, y_lower = read_data(
    filename=r'Opt_results_ideal',
    column=[2, 3, 4, 5]) * 1000000
no = len(wg_length)
# print wg_length

# Adjust the DC shape
pf.Exec("app.openproject(X:\BDC_pyfimm_August.10.prj,"")")
# It seems that the prj file have to be placed under the disk directly.

pf.Exec("app.subnodes[1].subnodes[1].lhsinput.inputtype = 2")
pf.Exec("app.subnodes[1].subnodes[1].lhsinput.setvec(0.5,0,0.5,0)")
pf.Exec("app.subnodes[1].subnodes[1].lhsinput.normalise = 1")


def taper_boundaries(w11, w12, w13, w31, w32, w33, l,
                     num):  # 2 additional variables need to be added: the lengths of input/output wgs
    z = np.linspace(0, l, 5)
    y1 = np.array(
        [WG_gap / 2. + WG_width, WG_gap / 2. + w11, WG_gap / 2. + w12, WG_gap / 2. + w13, WG_gap / 2. + WG_width])
    y3 = np.array(
        [-WG_gap / 2. - WG_width, -WG_gap / 2. - w31, -WG_gap / 2. - w32, -WG_gap / 2. - w33, -WG_gap / 2. - WG_width])

    z_new = np.linspace(z.min(), z.max(), num)
    y2 = np.full(num, WG_gap / 2.)
    y4 = np.full(num, -WG_gap / 2.)

    y1_cs = CubicSpline(z, y1, bc_type=((1, 0.0), (1, 0.0)))
    y3_cs = CubicSpline(z, y3, bc_type=((1, 0.0), (1, 0.0)))

    return z_new, y1_cs(z_new), y2, y3_cs(z_new), y4


for i in np.arange(0, no, 1):
    w11 = y_upper[i]
    w12 = y_upper[i]
    w13 = y_upper[i]
    w31 = y_lower[i]
    w32 = y_lower[i]
    w33 = y_lower[i]
    number = 50
    # print i
    # print wg_length[i]

    z_new, y1, y2, y3, y4 = taper_boundaries(w11, w12, w13, w31, w32, w33, T_length[i], number)

    for wl in np.linspace(1.5, 1.6, 11):
        pf.Exec("app.subnodes[1].subnodes[1].lambda=" + str(wl))
        pf.Exec("app.subnodes[1].subnodes[2].evlist.svp.lambda=" + str(wl))

        for j in range(number):
            w2 = y1[j] - y2[j]
            w4 = y4[j] - y3[j]
            w3 = WG_gap
            w1 = Clad_width + WG_width - w2
            w5 = Clad_width + WG_width - w4
            t_name = 'taper' + str(j)
            pf.Exec("Ref& {t_name}=app.subnodes[1].subnodes[" + str(2 + j) + "]")
            pf.Exec("app.subnodes[1].subnodes[" + str(j + 2) + "].evlist.svp.lambda=" + str(wl))
            pf.Exec("{t_name}.slices[1].width=" + str(w1))
            pf.Exec("{t_name}.slices[2].width=" + str(w2))
            pf.Exec("{t_name}.slices[3].width=" + str(w3))
            pf.Exec("{t_name}.slices[4].width=" + str(w4))
            pf.Exec("{t_name}.slices[5].width=" + str(w5))

        number = 50
        # DC_1
        pf.Exec("app.subnodes[1].subnodes[1].cdev.eltlist[1].length=" + str(wg_length[i]))
        # DC_2
        pf.Exec(
            "app.subnodes[1].subnodes[1].cdev.eltlist[" + str(2 * (number - 1) + 5) + "].length=" + str(wg_length[i]))

        # Taper/phase compensation section
        step_length = T_length[i] / number

        for k in range(number):
            pf.Exec("app.subnodes[1].subnodes[1].cdev.eltlist[" + str(2 * k + 3) + "].length=" + str(step_length))

        pf.Exec("app.subnodes[1].subnodes[1].update()")

        lr11 = pf.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[1][1]")
        lr12 = pf.Exec(
            "app.subnodes[1].subnodes[1].cdev.smat.lr[1][2]")  # ratio between amplitude of mode 2 at output and amplitude of mode 1 at input.
        lr21 = pf.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[2][1]")
        lr22 = pf.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[2][2]")

        Power_1 = 0.25 * np.abs((lr11 + lr21) + (lr12 + lr22)) ** 2
        Power_2 = 0.25 * np.abs((lr11 + lr21) - (lr12 + lr22)) ** 2

        Power_ratio = Power_1 / (Power_1 + Power_2)
        print "@", i, wl, wg_length[i], T_length[i], y_upper[i], y_lower[i], Power_ratio
