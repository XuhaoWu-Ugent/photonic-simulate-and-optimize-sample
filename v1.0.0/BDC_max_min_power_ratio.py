import pandas as pd
import numpy as np
import sys, os
from utils import *

#pyfimm = pd.read_csv("BDC_pyfimm_v1.csv")
#pyfimm.to_csv("BDC_pyfimm.csv", encoding='utf-8', index=False)

#for i in np.arange(2, 6, 1):
    #name = 'BDC_pyfimm' + str(i)
    #{name} = pd.read_csv("BDC_pyfimm"+ "_v" + str(i) + ".csv")
    #{name} = {name}[["wg_length", "T_length", "wavelength", "Power_ratio"]]
    #BDC_pyfimm = pd.read_csv("BDC_pyfimm.csv")
    #BDC_pyfimm = BDC_pyfimm[["wg_length", "T_length", "wavelength", "Power_ratio"]] hints for someone else for the loop in someday...
    #BDC_pyfimm = pd.concat([BDC_pyfimm, {name}], axis=0)
    
df1 = pd.read_csv("BDC_pyfimm_v1.csv")
df2 = pd.read_csv("BDC_pyfimm_v2.csv")
df3 = pd.read_csv("BDC_pyfimm_v3.csv")
df4 = pd.read_csv("BDC_pyfimm_v4.csv")
df5 = pd.read_csv("BDC_pyfimm_v5.csv")
df6 = pd.read_csv("BDC_pyfimm_v6.csv")
BDC_pyfimm = pd.concat([df1, df2, df3, df4, df5, df6], axis=0)
BDC_pyfimm = BDC_pyfimm[["wg_length", "T_length", "wavelength", "Power_ratio"]]
BDC_pyfimm.to_csv("BDC_pyfimm_total.csv", encoding='utf-8', index=False)

ver_num = 1
filename = "max_min_coupling.txt"
while os.path.isfile(filename):
    filename = "max_min_coupling"+"_v"+str(ver_num)+".txt"
    ver_num += 1
outfile = file(filename,'w')
print >> outfile, "wg_length T_length min_max_coupling average_power" 
    
for wg_length in np.linspace(10.0, 60.0, 101):
    
    BDC_pyfimm_wg = BDC_pyfimm[BDC_pyfimm['wg_length'] == wg_length]
    
    for T_length in np.linspace(5.0, 55.0, 101):
        BDC_pyfimm_t = BDC_pyfimm_wg[BDC_pyfimm_wg['T_length'] == T_length]
        
        Power_ratio_csv = BDC_pyfimm_t[["Power_ratio"]]
        Power_ratio_csv.to_csv("power_ratio_xuhao.csv", encoding='utf-8', index=False)
        Power_ratio = read_data(
            filename=r'power_ratio_xuhao',
            column=[0])        
        P_min = np.min(Power_ratio)
        P_max = np.max(Power_ratio)
        P_avg = np.average(Power_ratio)
        min_max_coupling = 10*np.log10(P_max/P_min) #also the average ratio (linear average, in dB)        
        print "@", wg_length, T_length, min_max_coupling, P_avg
        print >> outfile, wg_length, T_length, min_max_coupling, P_avg 
outfile.flush()
outfile.close()
        
print "finish"
        

