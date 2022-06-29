# Python module with simulation functions (cheap and expensive)\
# that can be called from Matlab in the sumo toolbox.

# This module is loaded by Matlab and runs in the same process.
# the functions call a simulation over RPC (in another process or on another PC)
# to execute the simulations.

print "Loading Python module"
import numpy as np
import sys
print sys.executable
import xmlrpclib
import time
from os import path


# Connect to the RPC server
s = xmlrpclib.ServerProxy('http://localhost:8000')

# initialize (will create the PCell)
s.siminit()

# variable to measure time (not CPU time)
end_time = time.time()

def cheap(x):
    global end_time
    t0 = time.time()
    t = s.DC_simulator([float(i) for i in x])
    delta_t = time.time()-t0
    if t>1:
        success = False
        t=1.0
    else:
        success = True

    add_to_file(simtype="CHEAP", result=t, variables=x, setup_time=t0-end_time, exec_time=delta_t, success=success)
    end_time = time.time()
    return float(t) # minimizer

def expensive(x):
    global end_time
        # check if there is existing data cached, not to rerun the same simulation
    ed = existing_data_generator.next()
    if ed is not None:
        print "Existing data"
        t = ed[1]
        success = ed[8]
        existing_x = ed[2:7]
        delta_t = ed[7]
        run = False #existing_x != x
    else:
        run = True

    if run:
        t0 = time.time()
        t = s.lumerical_fdtd_multi_wavelength([float(i) for i in x])
        delta_t = time.time()-t0
        if t>1:
            success = False
            t=1.0
        else:
            success = True

    add_to_file(simtype="EXPENSIVE", result=t, variables=x, setup_time=t0-end_time, exec_time=delta_t, success=success)
    end_time = time.time()
    return float(t) #minimizer

# file paths: file if needing to reload existing data (to continue simulation)
existing_file_path = r"test_results_yufei_existing.txt"
reload=False

# output file
file_path = r"test_results_xuhao.txt" #"test_results_yufei.txt"
# I am ashamed to change the file name, because I could not reproduce every step in this code list. Thanks to Yufei's great work, I could be here.

existing_data = []
if reload:
    # continue previous simulations: load data from file
    f = file(existing_file_path, "rt")
    lines = f.readlines()
    if len(lines)>1:
        for l in lines[1:]:
            if len(l) == 9:
                (typ, r, wpx, wpy, sbp, sbw, la, tim, success) = l.split(" ")
                sutim = "0.0"
            else:
                (typ, r, wpx, wpy, sbp, sbw, la, sutim, tim, success) = l.split(" ")
            existing_data.append((typ, float(r),
                                  float(wpx), float(wpy),
                                  float(sbp), float(sbw),
                                  float(la), float(sutim), float(tim),

                                  success=="True"))
            print "Added existing data line"
    f.close()

# wrap existing data in a generator object
def get_existing_data():
    for d in existing_data:
        yield d
    while 1:
        yield None

existing_data_generator = get_existing_data()

# open output file
f = file(file_path, "wt")
print >> f, "type result wg_length T_length y_upper y_lower setup_time sim_time success"

def add_to_file(simtype, result, variables, setup_time, exec_time, success):
    print >> f, simtype, result, 
    for v in variables:
        print >> f, v,
    print >> f, setup_time, exec_time, success
    f.flush()
    return 0


if __name__ == "__main__":
        # Testing the simulations
        # This is not executed in normal use from Matlab

    x = [30.2202,    0.0092]
    s.set_cell_params(x)
    s.visualize()
    print s.lumerical_fdtd_multi_wavelength(x, True)
    #print s.camfr_multi_wavelength(x, True, 100, "result_camfr.txt")
    #print s.cst_multi_wavelength(x, True, True)
    print s.get_times()
