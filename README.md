# photonic-simulate-and-optimize-sample
A photonic simulation and optimization sample based on my master thesis. The link code for other simulator and optimizer will be added later.



v 1.0.0

The supported simulator so far: Fimmwave

The supported optimizer so far: suno toolbox (Ugent)

How to use:

for BDC generation: BDC_pyfimm.

for parameter sweeping: BDC_sweep_length_pyfimm

for data analysis:
1. coupling ratio at a certain wl: contour_plot_pyfimm
2. max/min power: contour_plot_max_min_power
3. average power: contour_plot_average_power

for finding the potential candidates after sweeping: find_estimated_BDC

for Kriging running: pythonSimulator 

for finding the cadidates after the Kriging optimization: BDC_examination (for large number of points), BDC_examination_sample (for one point)
