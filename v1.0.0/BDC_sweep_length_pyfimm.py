import pyfimm as pf   # Every script must begin with this line
#pf.set_DEBUG()      # Enable Debugging output
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from time import gmtime, strftime, sleep
import pandas as pd

subprocess.Popen(r'"C:\Program Files (x86)\PhotonD6.6.2\Fimmwave\bin64\fimmwave.exe" -pt 5101') # Open fimmwave software.
sleep(10) # Wait for a couple of minutes to let App open.
print "Start"
pf.connect()        # this connects to the FimmWave application. 
begin_time = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
print "This program starts at {}".format(begin_time)

# Set Parameters (Your copy of FIMMWAVE has default values for these. You can change more than shown here. See __jaredwave.py
import sys, os
ScriptPath, ScriptFile = os.path.split( os.path.realpath(__file__)  )                    # Get directory of this script

pf.Exec("app.openproject(X:\BDC_pyfimm_June.6.prj,"")") 
# It seems that the prj file have to be placed under the disk directly.

# ############################################################
#sweep geometrical parameters: wavelength, waveguide length and taper length
#ver_num = 1
#filename = "BDC_pyfimm.txt"
#while os.path.isfile(filename):
    #filename = "BDC_pyfimm"+"_v"+str(ver_num)+".txt"
    #ver_num += 1
#outfile = file(filename,'w')
#print >> outfile, "wg_length T_length wl Power_1 Power_2 Power_ratio"

pf.Exec("app.subnodes[1].subnodes[1].lhsinput.inputtype = 2")
pf.Exec("app.subnodes[1].subnodes[1].lhsinput.setvec(0.5,0,0.5,0)")
pf.Exec("app.subnodes[1].subnodes[1].lhsinput.normalise = 1")

#adjust the shape of BDC

#widths
WG_width = 0.3
WG_gap = 0.2
Clad_width = 1
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

T_length = 5.0 # using float here to maintain the step length value is a float. (only necessary in python 2.7)
num_segments = 3 # The number of control points in taper section is 3*2=6.
w11 = 0.15
w12 = 0.15
w13 = 0.15
w31 = 0.45
w32 = 0.45
w33 = 0.45
num = 50
#strip_list = []
#T_strip = np.zero(num)

z_new, y1, y2, y3, y4 = taper_boundaries(w11, w12, w13, w31, w32, w33, T_length, num)   
plt.plot(z_new, y1, label='y1')
plt.plot(z_new, y2, label='y2')
plt.plot(z_new, y3, label='y3')
plt.plot(z_new, y4, label='y4')
plt.legend()
plt.show() # the derivative of curves at start/end points = 0.
# remember to close the fig window afterwards to let the code continue to run.

for i in range(num):
    w2 = y1[i] - y2[i]
    w4 = y4[i] - y3[i]
    w3 = WG_gap
    w1 = Clad_width + WG_width - w2
    w5 = Clad_width + WG_width - w4
    t_name = 'taper' + str(i)
    pf.Exec("Ref& {t_name}=app.subnodes[1].subnodes["+str(2+i)+"]")
    pf.Exec("{t_name}.slices[1].width="+str(w1))
    pf.Exec("{t_name}.slices[2].width="+str(w2))
    pf.Exec("{t_name}.slices[3].width="+str(w3))
    pf.Exec("{t_name}.slices[4].width="+str(w4))
    pf.Exec("{t_name}.slices[5].width="+str(w5))

for wl in np.linspace(1.5,1.6,6):
    ver_num = 1
    filename = "BDC_pyfimm.txt"
    while os.path.isfile(filename):
        filename = "BDC_pyfimm"+"_v"+str(ver_num)+".txt"
        ver_num += 1
    outfile = file(filename,'w')
    print >> outfile, "wg_length T_length wavelength Power_1 Power_2 Power_ratio"
    
    pf.Exec("app.subnodes[1].subnodes[1].lambda="+str(wl))     
    
    for wg_length in np.linspace(8.0, 20.0, 41): #more points
        num = 50
        pf.Exec("app.subnodes[1].subnodes[2].evlist.svp.lambda="+str(wl))
        # DC_1
        pf.Exec("app.subnodes[1].subnodes[1].cdev.eltlist[1].length="+str(wg_length))
        # DC_2
        pf.Exec("app.subnodes[1].subnodes[1].cdev.eltlist["+str(2*(num-1)+5)+"].length="+str(wg_length))  
        
        for T_length in np.linspace(16.0, 55.0, 41): #more points, longer length
            step_length = step_length = T_length/num
            
            for i in range(num):
                pf.Exec("app.subnodes[1].subnodes["+str(i+3)+"].evlist.svp.lambda="+str(wl))
                pf.Exec("app.subnodes[1].subnodes[1].cdev.eltlist["+str(2*i+3)+"].length="+str(step_length))
                
            pf.Exec("app.subnodes[1].subnodes[1].update()") #note the indent
            lr11 = pf.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[1][1]") 
            lr12 = pf.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[1][2]") #ratio between amplitude of mode 2 at output and amplitude of mode 1 at input.
            lr21 = pf.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[2][1]")
            lr22 = pf.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[2][2]")
            #Power_1 = 0.5*(lr11*np.conjugate(lr11) + lr21*np.conjugate(lr21))
            Power_1 = 0.25*np.abs((lr11 + lr21) + (lr12 + lr22))**2
            Power_2 = 0.25*np.abs((lr11 + lr21) - (lr12 + lr22))**2 #doing the calculation, check the sign
            #Power_2 = 0.5*(lr12*np.conjugate(lr12) + lr22*np.conjugate(lr22))
            Power_ratio = Power_1/(Power_1+Power_2) #
                
            print "@", wl, wg_length, T_length, lr11, lr12, lr21, lr22, Power_1, Power_2, Power_ratio
            print >> outfile, wg_length, T_length, wl, Power_1, Power_2, Power_ratio
    outfile.flush()
    outfile.close()
 
    df = pd.read_csv( "BDC_pyfimm"+"_v"+str(ver_num-1)+".txt",delimiter=" ")
    df.to_csv( "BDC_pyfimm"+"_v"+str(ver_num-1)+".csv", encoding='utf-8', index=False)

