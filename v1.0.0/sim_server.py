# Very Simulation Server
# Launch this script as a separate process. 
# as long as it runs it will operate an RPC server
# that launches simulations when triggered over the network

# This allows you to run simulations on different PCs,
# or in different processes on the same PC 
# e.g. Matlab and CAMFR cannot run in the same process.

from __future__ import division
import time
# simple simulation server
# every method can be called over RPC
# just add different methods for each simulation you would like to run
# (Hey, no-one said I have to make this generic.)
from utils import *
import subprocess
import numpy as np
from pdPythonLib import *
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline


class SimServer(object):

    def siminit(self):
        """ Initialize:
            Load IPKISS and the cell to be simulated
        """
        print "Loading Python module"
        import sys
        print sys.executable

        self.reset()
        return 0

    def reset(self):
        """ reset all execution statistics """
        
        self.cheap_times = []
        self.expensive_times = []
        self.cheap_fail_times = []
        self.expensive_fail_times = []
        self.cheap_fails = []
        self.expensive_fails = []

    def set_cell_params(self,x):
        """ Set the parameters of the simulated Cell.
            
            x: np.array
            x[0] coupler length
        """ 

        return 0
    
    def dummy_cheap(self, x, show=False):
        import numpy as np
        return float(0.5*np.mean(x))

    def get_times(self):
        return [self.cheap_times, self.expensive_times, self.cheap_fail_times, self.expensive_fail_times]

    def get_failures(self):
        return [self.cheap_fails, self.expensive_fails]    
            
    def DC_simulator(self,x):
        print x
        p = subprocess.Popen(r'"C:\Program Files (x86)\PhotonD6.6.2\Fimmwave\bin64\fimmwave.exe" -pt 5101') # Open fimmwave software.
        time.sleep(10) 
        print "start"
        
        fimm = pdApp()
        fimm.ConnectToApp()                      
        
        fimm.Exec("app.openproject(X:\BDC_pyfimm_June.6.prj,"")") 
        
        fimm.Exec("app.subnodes[1].subnodes[1].lhsinput.inputtype = 2")
        fimm.Exec("app.subnodes[1].subnodes[1].lhsinput.setvec(0.5,0,0.5,0)")
        fimm.Exec("app.subnodes[1].subnodes[1].lhsinput.normalise = 1")  
        
        #widths
        WG_width = 0.3
        WG_gap = 0.2
        Clad_width = 1        
        
        num_segments = 3 # The number of control points in taper section is 3*2=6.       
        
        num_wl = 5
        wg_length = x[0]*1000000
        T_length = x[1]*1000000
        w11 = x[2]*1000000
        w12 = x[2]*1000000
        w13 = x[2]*1000000
        w31 = x[3]*1000000
        w32 = x[3]*1000000
        w33 = x[3]*1000000
        
        
        print wg_length, T_length, w11, w31
        P_out = []
        
        def taper_boundaries(w11, w12, w13, w31, w32, w33, l, num): #this part of code is for defining the shape taper/PC section in DC_simulator.
            z = np.linspace(0, l, 5)
            y1 = np.array([WG_gap/2.+WG_width,WG_gap/2.+w11, WG_gap/2.+w12, WG_gap/2.+w13, WG_gap/2.+WG_width])
            y3 = np.array([-WG_gap/2.-WG_width,-WG_gap/2.-w31, -WG_gap/2.-w32, -WG_gap/2.-w33, -WG_gap/2.-WG_width])
        
            z_new = np.linspace(z.min(), z.max(),num)
            y2 = np.full(num, WG_gap/2.)
            y4 = np.full(num, -WG_gap/2.)
        
            y1_cs = CubicSpline(z, y1, bc_type=((1, 0.0), (1, 0.0)))
            y3_cs = CubicSpline(z, y3, bc_type=((1, 0.0), (1, 0.0)))
            
            return z_new, y1_cs(z_new), y2, y3_cs(z_new), y4        
        
        num = 50
        
        z_new, y1, y2, y3, y4 = taper_boundaries(w11, w12, w13, w31, w32, w33, T_length, num) 
        
        for i in range(num):
            w2 = y1[i] - y2[i]
            w4 = y4[i] - y3[i]
            w3 = WG_gap
            w1 = Clad_width + WG_width - w2
            w5 = Clad_width + WG_width - w4
            t_name = 'taper' + str(i)
            fimm.Exec("Ref& {t_name}=app.subnodes[1].subnodes[".format(t_name=t_name)+str(2+i)+"]") #note that the variable "t_name" need to be formated, otherwise it will be regarded as undefined.
            fimm.Exec("{t_name}.slices[1].width=".format(t_name=t_name)+str(w1))
            fimm.Exec("{t_name}.slices[2].width=".format(t_name=t_name)+str(w2))
            fimm.Exec("{t_name}.slices[3].width=".format(t_name=t_name)+str(w3))
            fimm.Exec("{t_name}.slices[4].width=".format(t_name=t_name)+str(w4))
            fimm.Exec("{t_name}.slices[5].width=".format(t_name=t_name)+str(w5))        
        
        for wl in np.linspace(1.5,1.58,num_wl): 
            fimm.Exec("app.subnodes[1].subnodes[1].lambda="+str(wl))
            
            num = 50
            fimm.Exec("app.subnodes[1].subnodes[2].evlist.svp.lambda="+str(wl)) #change the wavelength in cross section of waveguide
            fimm.Exec("app.subnodes[1].subnodes[1].cdev.eltlist[1].length="+str(wg_length))
            fimm.Exec("app.subnodes[1].subnodes[1].cdev.eltlist["+str(2*(num-1)+5)+"].length="+str(wg_length)) 
            step_length = T_length/num
            
            for i in range(num):
                fimm.Exec("app.subnodes[1].subnodes["+str(i+3)+"].evlist.svp.lambda="+str(wl)) #change the wavelength in every cross section of taper/PC
                fimm.Exec("app.subnodes[1].subnodes[1].cdev.eltlist["+str(2*i+3)+"].length="+str(step_length))
                
            fimm.Exec("app.subnodes[1].subnodes[1].update()") #note the indent
            lr11 = fimm.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[1][1]")
            lr12 = fimm.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[1][2]")
            lr21 = fimm.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[2][1]")
            lr22 = fimm.Exec("app.subnodes[1].subnodes[1].cdev.smat.lr[2][2]")
            Power_1 = 0.25*np.abs((lr11 + lr21) + (lr12 + lr22))**2
            Power_2 = 0.25*np.abs((lr11 + lr21) - (lr12 + lr22))**2
            Power_ratio = Power_1/(Power_1 + Power_2)
            print Power_ratio
            P_out.append(Power_ratio)
            
        p.terminate()
        return float(max(abs(P_out-np.array([0.50]*num_wl)))) # change the coupling here 
        

if __name__ == "__main__":
    # Launch the server
    from SimpleXMLRPCServer import SimpleXMLRPCServer
    from SimpleXMLRPCServer import SimpleXMLRPCRequestHandler    
    # Restrict to a particular path.
    class RequestHandler(SimpleXMLRPCRequestHandler):
        rpc_paths = ('/RPC2',)    

    # Create server
    server = SimpleXMLRPCServer(("localhost", 8000),
                                requestHandler=RequestHandler)
    server.register_introspection_functions()

    server.register_instance(SimServer())

    server.serve_forever()



    
