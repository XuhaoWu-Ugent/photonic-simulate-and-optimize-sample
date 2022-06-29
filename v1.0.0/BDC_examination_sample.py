import pyfimm as pf   # Every script must begin with this line
#pf.set_DEBUG()      # Enable Debugging output
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from time import gmtime, strftime, sleep
import pandas as pd
from utils import *

subprocess.Popen(r'"C:\Program Files (x86)\PhotonD6.6.2\Fimmwave\bin64\fimmwave.exe" -pt 5101') # Open fimmwave software.
sleep(10) # Wait for a couple of minutes to let App open.
print "Start"
pf.connect()        # this connects to the FimmWave application. 
begin_time = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
print "This program starts at {}".format(begin_time)

# Set Parameters (Your copy of FIMMWAVE has default values for these. You can change more than shown here. See __jaredwave.py
import sys, os
ScriptPath, ScriptFile = os.path.split( os.path.realpath(__file__)  )                    # Get directory of this script

#widths
WG_width = 0.3
WG_gap = 0.2
Clad_width = 1

#Adjust the DC shape
pf.Exec("app.openproject(X:\BDC_pyfimm_June.6.prj,"")") 
# It seems that the prj file have to be placed under the disk directly.

pf.Exec("app.subnodes[1].subnodes[1].lhsinput.inputtype = 2")
pf.Exec("app.subnodes[1].subnodes[1].lhsinput.setvec(0.5,0,0.5,0)")
pf.Exec("app.subnodes[1].subnodes[1].lhsinput.normalise = 1")

def taper_boundaries(w11, w12, w13, w31, w32, w33, l, num): #2 additional variables need to be added: the lengths of input/output wgs
    z = np.linspace(0, l, 5)
    y1 = np.array([WG_gap/2.+WG_width,WG_gap/2.+w11, WG_gap/2.+w12, WG_gap/2.+w13, WG_gap/2.+WG_width])
    y3 = np.array([-WG_gap/2.-WG_width,-WG_gap/2.-w31, -WG_gap/2.-w32, -WG_gap/2.-w33, -WG_gap/2.-WG_width])

    z_new = np.linspace(z.min(), z.max(),num)
    y2 = np.full(num, WG_gap/2.)
    y4 = np.full(num, -WG_gap/2.)

    y1_cs = CubicSpline(z, y1, bc_type=((1, 0.0), (1, 0.0)))
    y3_cs = CubicSpline(z, y3, bc_type=((1, 0.0), (1, 0.0)))
    
    return z_new, y1_cs(z_new), y2, y3_cs(z_new), y4
#7.336740057889999e-06,1.55107489662e-05,2.1572706357600003e-07,5.487910909270001e-07
wg_length = 7.336740057889999 #the parameters of one sample point
T_length = 15.5107489662
w11 = 0.21572706357600003
w12 = 0.21572706357600003
w13 = 0.21572706357600003
w31 = 0.5487910909270001
w32 = 0.5487910909270001
w33 = 0.5487910909270001
number = 50
power_ratio = []
num_wl = 81

z_new, y1, y2, y3, y4 = taper_boundaries(w11, w12, w13, w31, w32, w33, T_length, number)       

for wl in np.linspace(1.5,1.58,num_wl):
    pf.Exec("app.subnodes[1].subnodes[1].lambda="+str(wl))
    pf.Exec("app.subnodes[1].subnodes[2].evlist.svp.lambda="+str(wl))
    
    for i in range(number):
        w2 = y1[i] - y2[i]
        w4 = y4[i] - y3[i]
        w3 = WG_gap
        w1 = Clad_width + WG_width - w2
        w5 = Clad_width + WG_width - w4
        t_name = 'taper' + str(i)
        pf.Exec("Ref& {t_name}=app.subnodes[1].subnodes["+str(2+i)+"]")
        pf.Exec("app.subnodes[1].subnodes["+str(i+2)+"].evlist.svp.lambda="+str(wl))
        pf.Exec("{t_name}.slices[1].width="+str(w1))
        pf.Exec("{t_name}.slices[2].width="+str(w2))
        pf.Exec("{t_name}.slices[3].width="+str(w3))
        pf.Exec("{t_name}.slices[4].width="+str(w4))
        pf.Exec("{t_name}.slices[5].width="+str(w5)) 
        
    number = 50
    # DC_1
    pf.Exec("app.subnodes[1].subnodes[1].cdev.eltlist[1].length="+str(wg_length))
    # DC_2
    pf.Exec("app.subnodes[1].subnodes[1].cdev.eltlist["+str(2*(number-1)+5)+"].length="+str(wg_length))  
            
    #Taper/phase compensation section
    step_length = T_length/number
                
    for j in range(number):
        pf.Exec("app.subnodes[1].subnodes[1].cdev.eltlist["+str(2*j+3)+"].length="+str(step_length))   
        
    pf.Exec("app.subnodes[1].subnodes[1].update()") 
        
    lr11 = pf.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[1][1]") 
    lr12 = pf.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[1][2]") #ratio between amplitude of mode 2 at output and amplitude of mode 1 at input.
    lr21 = pf.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[2][1]")
    lr22 = pf.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[2][2]")
        
    Power_1 = 0.25*np.abs((lr11 + lr21) + (lr12 + lr22))**2
    Power_2 = 0.25*np.abs((lr11 + lr21) - (lr12 + lr22))**2 
        
    Power_ratio = Power_1/(Power_1+Power_2) 
    print "@", wl, Power_ratio
    power_ratio.append(Power_ratio)
    
# make the plot
wavelength = np.linspace(1.5,1.58,num_wl)
power_ratio = np.array(power_ratio)
wavelength = np.array(wavelength) # these "np.array" are necessary, otherwise ValueError: x and y must have same first dimension will occurs. a bug of np acutually.
P_avg = np.average(power_ratio)
P_max_min = 10*np.log10(np.max(power_ratio)/np.min(power_ratio))
plt.plot(wavelength*1000, power_ratio*100, "b-o")
plt.xlabel(r'wavelength [nm]')
plt.ylabel(r'Coupling power ratio [%]')
plt.title(r'The coupling power ratio for $lambda$ in 1500-1580 nm')
print "average_power = ", P_avg, "ratio_of_max_to_min_power = ", P_max_min
plt.show()